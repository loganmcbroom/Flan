#include "flan/Wavetable.h"

#include <algorithm>
#include <numeric>
#include <iostream>
//#include <ranges>

#include "WDL/resample.h"
#include "r8brain/CDSPResampler.h"
#include "flan/PV/PV.h"
#include "flan/Graph.h"

#undef min
#undef max

using namespace flan;

static Frame snap_frame_to_sample( const Audio & me, Frame frameToSnap, Sample snapHeight, Frame searchDistance )
	{	
	const Frame leftFrameLimit = std::max( frameToSnap - searchDistance, 0 );
	const Frame rightFrameLimit = std::min( frameToSnap + searchDistance, me.get_num_frames() - 1 );

	// Bidirectionally search from start_frame until a crossing is found
	const bool isAbove = me.get_sample( 0, frameToSnap ) > snapHeight;
	for( Frame searchOffset = 0; searchOffset <= searchDistance; ++searchOffset )
		{
		// Check left
		const Frame leftSearchFrame = frameToSnap - searchOffset;
		if( leftSearchFrame >= leftFrameLimit && me.get_sample( 0, leftSearchFrame ) > snapHeight != isAbove )
			return leftSearchFrame + 1;
			
		// Check right
		const Frame rightSearchFrame = frameToSnap + searchOffset;
		if( rightSearchFrame < rightFrameLimit && me.get_sample( 0, rightSearchFrame ) > snapHeight != isAbove )
			return rightSearchFrame;
		}
	
	// If control reaches here, cross search failed. Use frame with output nearest to crossing.
	auto norm = [snapHeight, &me, frameToSnap, searchDistance]( Frame f )
		{ 
		const float m = 1.0f;
		const float r = 1.0f + m * float( std::abs( f - frameToSnap ) ) / searchDistance; // Slightly prefer frames nearer to initial
		return std::abs( me.get_sample( 0, f ) - snapHeight ) * r;
		};

	Frame currentBestFrame = frameToSnap;
	Sample currentBestSampleDistance = norm( frameToSnap );
	for( Frame searchFrame = leftFrameLimit; searchFrame <= rightFrameLimit; ++searchFrame )
		{
		const Sample currentSampleDistance = norm( searchFrame );
		if( currentSampleDistance < currentBestSampleDistance )
			{
			currentBestFrame = searchFrame;
			currentBestSampleDistance = currentSampleDistance;
			}
		}

	return currentBestFrame;
	}

static Audio resample_waveforms( const Audio & source, const std::vector<Frame> & waveform_starts, Frame outputWavelength, flan_CANCEL_ARG_CPP )
    {
	// Input validation
	if( source.is_null() ) return Audio::create_null();
	if( waveform_starts.empty() ) return Audio::create_null();
	const Audio mono_source = source.get_num_channels() == 1? source.copy() : source.convert_to_mono();

	// Set up output
    Audio::Format format = mono_source.get_format();
	format.num_channels = 1;
    format.num_frames = outputWavelength * ( waveform_starts.size() - 1 );
	Audio out( format );

	// For each waveform, resample to wavelength
	for( Waveform waveform = 0; waveform < waveform_starts.size() - 1; ++waveform )
		{
		flan_CANCEL_POINT( Audio::create_null() );

		const Frame start_frame = waveform_starts[waveform];
		const Frame end_frame =  waveform_starts[waveform + 1];
		const Frame numInputFrames = end_frame - start_frame;
		const float expansion = float( outputWavelength ) / numInputFrames;
		const Frame paddedStartFrame = std::max( start_frame - numInputFrames, 0 );
		const Frame paddedEndFrame = std::min( end_frame + numInputFrames, mono_source.get_num_frames() );
		const Frame paddedLength = paddedEndFrame - paddedStartFrame;
		const Frame leftPadSize = start_frame - paddedStartFrame;
		const Frame rightPadSize = paddedEndFrame - end_frame;

		// Resample
		r8b::CDSPResampler rs( numInputFrames, outputWavelength, 3 * numInputFrames );
		const Sample * inPtr = mono_source.get_sample_pointer( 0, paddedStartFrame );
		std::vector<Sample> outBuffer( paddedLength * expansion );
		rs.oneshot( inPtr, paddedLength, outBuffer.data(), outBuffer.size() );
	
		std::copy( FLAN_PAR_UNSEQ
			outBuffer.begin() + leftPadSize * expansion, 
			outBuffer.begin() + leftPadSize * expansion + outputWavelength, 
			out.get_sample_pointer( 0, waveform * outputWavelength ) );
		}

    return out;
    }

std::vector<Frame> Wavetable::get_waveform_starts( 
	const Audio & source, 
	Wavetable::SnapMode snap_mode, 
	Wavetable::PitchMode pitch_mode, 
	float snap_ratio, 
	Frame fixed_frame, 
	flan_CANCEL_ARG_CPP )
	{
	// Input validation
	if( source.is_null() ) return std::vector<Frame>();
	const Audio mono_source = source.convert_to_mono(); // Convert input to mono

	// Find expected waveform lengths via autocorrelation-based pitch detection
	const Frame ac_granularity = 128;
	const Frame ac_windowSize = 4096;
	std::vector<fFrame> local_wavelengths;
	Frame global_wavelength = 0;
	if( pitch_mode != Wavetable::PitchMode::None )
		{
		local_wavelengths = mono_source.get_local_wavelengths( 0, 0, -1, 2048, 128, canceller );
		global_wavelength = mono_source.get_average_wavelength( local_wavelengths, .5, 64 ); 
		if( global_wavelength == -1 ) 
			{
			std::cout << "Insignificant pitch data found in Wavetable construction.\n";
			return std::vector<Frame>();
			}
		}

	std::vector<Frame> waveform_starts;

	auto snap_handler = [&]( Frame frameToSnap, Frame snapSourceFrame, Frame maxSnap )
		{
		//const Frame maxSnap = snap_ratio * std::abs( frameToSnap - snapSourceFrame );
		switch( snap_mode )
			{
			case Wavetable::SnapMode::None: return frameToSnap;
			case Wavetable::SnapMode::Zero: return snap_frame_to_sample( mono_source, frameToSnap, 0, maxSnap );
			case Wavetable::SnapMode::Level: return snap_frame_to_sample( mono_source, frameToSnap, mono_source.get_sample( 0, snapSourceFrame ), maxSnap );
			default: return frameToSnap;
			}
		};
	
	waveform_starts.push_back( snap_handler( 0, 0, snap_ratio * global_wavelength ) );

	while( true ) // Eat until we run out of buffer
		{
		flan_CANCEL_POINT( std::vector<Frame>() );

		// Guess how many frames the next wavelength will be
		Frame expected_num_frames = 0;
		if( pitch_mode ==  Wavetable::PitchMode::Local )
			{
			// Where in the local wavelegnths vector we currently are
			const int local_wavelength_index = std::floor( waveform_starts.back() / ac_granularity );
			if( local_wavelength_index >= local_wavelengths.size() ) break;
			const Frame local_wavelength_c = local_wavelengths[local_wavelength_index];
			expected_num_frames = local_wavelength_c > 0? local_wavelength_c : global_wavelength; 
			}
		else if( pitch_mode == Wavetable::PitchMode::Global )
				expected_num_frames = global_wavelength; 
		else if( pitch_mode == Wavetable::PitchMode::None )
			expected_num_frames = fixed_frame;

		if( waveform_starts.back() + expected_num_frames >= mono_source.get_num_frames() ) break;

		// Snap expected waveform end if needed
		waveform_starts.push_back( snap_handler( waveform_starts.back() + expected_num_frames, waveform_starts.back(), snap_ratio * expected_num_frames ) );
		}

    return waveform_starts;
	}

Wavetable::Wavetable( const Audio & source, SnapMode snap_mode, PitchMode pitch_mode, float snap_ratio, Frame fixed, flan_CANCEL_ARG_CPP )
	: wavelength( 2048 )
    , table( resample_waveforms( source, get_waveform_starts( source, snap_mode, pitch_mode, snap_ratio, fixed, canceller ), wavelength, canceller ).set_volume( 1 ) )
	, num_waveforms( table.get_num_frames() / wavelength )
	{
	}

Wavetable::Wavetable( const Function<Second, Amplitude> & f, int num_waves, flan_CANCEL_ARG_CPP )
	: wavelength( 2048 )
    , table()
	, num_waveforms( num_waves )
	{
	const int oversample = 16;

	// Set up output
    Audio::Format format;
	format.num_channels = 1;
    format.num_frames = wavelength * num_waves;
	table = Audio( format );

	const int overLength = wavelength * oversample;
	std::vector<Sample> fBuffer( overLength * 3 );

	r8b::CDSPResampler rs( overLength, wavelength, fBuffer.size() );

	// For each waveform, resample to wavelength
	for( Waveform waveform = 0; waveform < num_waves; ++waveform )
		{
		if( canceller ) return;
		// Sample f into buffer
		for( int s = 0; s < overLength; ++s )
			{
			fBuffer[s] = f( waveform + float( s ) / overLength );
			fBuffer[s + 1 * wavelength * oversample] = fBuffer[s];
			fBuffer[s + 2 * wavelength * oversample] = fBuffer[s];
			}

		// Resample
		rs.oneshot( fBuffer.data(), fBuffer.size(), table.get_sample_pointer( 0, wavelength * waveform ), wavelength );
		}
	}

Wavetable::Wavetable( const Audio & source, const std::vector<Frame> & waveform_starts, flan_CANCEL_ARG_CPP )
	: wavelength( 2048 )
    , table( resample_waveforms( source, waveform_starts, wavelength, canceller ).set_volume( 1 ) )
	, num_waveforms( waveform_starts.size() - 1 )
	{
	}

Wavetable::Wavetable()
	: wavelength( 0 )
    , table( Audio::create_null() )
	, num_waveforms( 0 )
	{
	}

Audio Wavetable::synthesize( 
	Second length, 
	const Function<Second, Frequency> & freq, 
	const Function<Second, float> & index, 
	Second granularity_time, 
	flan_CANCEL_ARG_CPP 
	) const
    {
	if( num_waveforms == 0 ) return Audio::create_null();
	
	// Ouput setup
	Audio::Format format = table.get_format();
    format.num_frames = table.time_to_frame( length );
    Audio out( format );

	const Frame granularity = out.time_to_frame( granularity_time );

	// Resampler setup
	WDL_Resampler rs;
	rs.SetMode( true, 1, true, 512 );
	
	Frame out_frames_generated = 0;
	Frame phaseFrame = 0;
	while( out_frames_generated < out.get_num_frames() )
		{
		flan_CANCEL_POINT( Audio::create_null() );

		const double inFreq_c = double( table.get_sample_rate() ) / wavelength;
		const double outFreq_c = freq( out.frame_to_time( out_frames_generated ) );
		const float index_c = std::clamp( num_waveforms * index( out.frame_to_time( out_frames_generated ) ), 0.0f, float( num_waveforms - 1 ) );
		const Frame leftIndex = std::floor( index_c );
		const Frame rightIndex = std::ceil( index_c );
		const float indexRemainder = index_c - leftIndex;

		// Current resample setup
		rs.SetRates( outFreq_c, inFreq_c );
		WDL_ResampleSample * rsinbuf = nullptr;
		const Frame wdlWanted = rs.ResamplePrepare( granularity, 1, &rsinbuf );

		// Copy waveform to wdl
		for( Frame j = 0; j < wdlWanted; ++j )
			{
			const Sample leftSample  = table.get_sample( 0, ( phaseFrame + j ) % wavelength + wavelength * leftIndex  );
			const Sample rightsample = table.get_sample( 0, ( phaseFrame + j ) % wavelength + wavelength * rightIndex );
			rsinbuf[j] = ( 1.0f - indexRemainder ) * leftSample + indexRemainder * rightsample;
			}
				
		// Resample directly into output buffer and update frame indices.
		out_frames_generated += rs.ResampleOut( out.get_sample_pointer( 0, out_frames_generated ), wdlWanted, out.get_num_frames() - out_frames_generated, 1 );
		phaseFrame = ( phaseFrame + wdlWanted ) % wavelength;
		}
	
	return out;
    }

void Wavetable::add_fades_in_place( Frame fade_frames, Waveform waveform )
	{
	const int waveform_index = std::clamp( int( waveform * num_waveforms ), 0, num_waveforms - 1 );

	auto singleTransform = [this, fade_frames]( int w )
		{
		Sample * waveform = get_waveform( w );
		for( Frame frame = 0; frame < fade_frames - 1; ++frame )
			{
			const float fade = std::sin( pi / 2.0f * ( frame + 1 ) / fade_frames );
			*(waveform + frame) *= fade;
			*(waveform + wavelength - 1 - frame) *= fade;
			}
		};

	if( waveform == -1 ) 
		for( int i = 0; i < num_waveforms; ++i )
			singleTransform( i );
	else singleTransform( waveform_index );
	}

void Wavetable::remove_jumps_in_place( Frame fade_frames, Waveform waveform )
	{
	const int waveform_index = std::clamp( int( waveform * num_waveforms ), 0, num_waveforms - 1 );

	auto singleTransform = [this, fade_frames]( int w )
		{
		Sample * waveform = get_waveform( w );
		const Sample l = *waveform;
		const Sample r = *( waveform + wavelength - 1 );
		const Sample mid = ( l + r ) / 2.0f;

		auto fader = [mid]( Sample * s, float fade ){ *s = ( *s - mid ) * fade + mid; };

		for( Frame frame = 0; frame < fade_frames - 1; ++frame )
			{
			const float fade = std::sin( pi / 2.0f * ( frame + 1 ) / fade_frames );
			fader( waveform + frame, fade );
			fader( waveform + wavelength - 1 - frame, fade );
			}
		};

	if( waveform == -1 ) 
		for( int i = 0; i < num_waveforms; ++i )
			singleTransform( i );
	else singleTransform( waveform_index );
	}

void Wavetable::remove_dc_in_place( Waveform waveform )
	{
	const int waveform_index = std::clamp( int( waveform * num_waveforms ), 0, num_waveforms - 1 );

	auto singleTransform = [this]( int w )
		{
		auto start_iter = table.get_buffer().begin() + wavelength * w;
		const float dc = std::accumulate( start_iter, start_iter + wavelength, 0.0f ) / float( wavelength );
		std::for_each( FLAN_PAR_UNSEQ start_iter, start_iter + wavelength, [dc]( Sample & s ){ s -= dc; } );
		};

	if( waveform == -1 ) 
		for( int i = 0; i < num_waveforms; ++i )
			singleTransform( i );
	else singleTransform( waveform_index );
	}

Sample * Wavetable::get_waveform( Waveform waveform )
	{
	const int waveform_index = std::clamp( int( waveform * num_waveforms ), 0, num_waveforms - 1 );
	return table.get_sample_pointer( 0, waveform_index * wavelength );
	}
	