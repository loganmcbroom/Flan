#include "flan/Wavetable.h"

#include <algorithm>
#include <numeric>
#include <iostream>
//#include <ranges>

#include "WDL/resample.h"
#include "r8brain/CDSPResampler.h"
#include "flan/PVOC.h"
#include "flan/Graph.h"

#undef min
#undef max

using namespace flan;

static const float pi = std::acos( -1.0f );
static const float pi2 = 2.0f * pi;

static Frame snapFrameToSample( const Audio & me, Frame frameToSnap, Sample snapHeight, Frame searchDistance, flan_CANCEL_ARG_CPP )
	{	
	const Frame leftFrameLimit = std::max( frameToSnap - searchDistance, 0 );
	const Frame rightFrameLimit = std::min( frameToSnap + searchDistance, me.getNumFrames() - 1 );

	// Bidirectionally search from startFrame until a crossing is found
	const bool isAbove = me.getSample( 0, frameToSnap ) > snapHeight;
	for( Frame searchOffset = 0; searchOffset <= searchDistance; ++searchOffset )
		{
		flan_CANCEL_POINT( 0 );

		// Check left
		const Frame leftSearchFrame = frameToSnap - searchOffset;
		if( leftSearchFrame >= leftFrameLimit && me.getSample( 0, leftSearchFrame ) > snapHeight != isAbove )
			return leftSearchFrame + 1;
			
		// Check right
		const Frame rightSearchFrame = frameToSnap + searchOffset;
		if( rightSearchFrame < rightFrameLimit && me.getSample( 0, rightSearchFrame ) > snapHeight != isAbove )
			return rightSearchFrame;
		}
	
	// If control reaches here, cross search failed. Use frame with output nearest to crossing.
	auto norm = [snapHeight, &me, frameToSnap, searchDistance]( Frame f )
		{ 
		const float m = 1.0f;
		const float r = 1.0f + m * float( std::abs( f - frameToSnap ) ) / searchDistance; // Slightly prefer frames nearer to initial
		return std::abs( me.getSample( 0, f ) - snapHeight ) * r;
		};

	Frame currentBestFrame = frameToSnap;
	Sample currentBestSampleDistance = norm( frameToSnap );
	for( Frame searchFrame = leftFrameLimit; searchFrame <= rightFrameLimit; ++searchFrame )
		{
		flan_CANCEL_POINT( 0 );

		const Sample currentSampleDistance = norm( searchFrame );
		if( currentSampleDistance < currentBestSampleDistance )
			{
			currentBestFrame = searchFrame;
			currentBestSampleDistance = currentSampleDistance;
			}
		}

	return currentBestFrame;
	}

static Audio resampleWaveforms( const Audio & source, const std::vector<Frame> & waveformStarts, Frame outputWavelength, flan_CANCEL_ARG_CPP )
    {
	// Input validation
	if( source.isNull() ) return Audio();
	const Audio monoSource = source.getNumChannels() == 1? source : source.convertToMono();

	// Set up output
    Audio::Format format = monoSource.getFormat();
	format.numChannels = 1;
    format.numFrames = outputWavelength * ( waveformStarts.size() - 1 );
	Audio out( format );

	// For each waveform, resample to wavelength
	for( Waveform waveform = 0; waveform < waveformStarts.size() - 1; ++waveform )
		{
		flan_CANCEL_POINT( Audio() );

		const Frame startFrame = waveformStarts[waveform];
		const Frame endFrame =  waveformStarts[waveform + 1];
		const Frame numInputFrames = endFrame - startFrame;
		const float expansion = float( outputWavelength ) / numInputFrames;
		const Frame paddedStartFrame = std::max( startFrame - numInputFrames, 0 );
		const Frame paddedEndFrame = std::min( endFrame + numInputFrames, monoSource.getNumFrames() );
		const Frame paddedLength = paddedEndFrame - paddedStartFrame;
		const Frame leftPadSize = startFrame - paddedStartFrame;
		const Frame rightPadSize = paddedEndFrame - endFrame;

		// Resample
		r8b::CDSPResampler rs( numInputFrames, outputWavelength, 3 * numInputFrames );
		const Sample * inPtr = monoSource.getSamplePointer( 0, paddedStartFrame );
		std::vector<Sample> outBuffer( paddedLength * expansion );
		rs.oneshot( inPtr, paddedLength, outBuffer.data(), outBuffer.size() );
	
		std::copy( outBuffer.begin() + leftPadSize * expansion, outBuffer.begin() + leftPadSize * expansion + outputWavelength, out.getSamplePointer( 0, waveform * outputWavelength ) );
		}

    return out;
    }

std::vector<Frame> Wavetable::getWaveformStarts( const Audio & source, Wavetable::SnapMode snapMode, 
	Wavetable::PitchMode pitchMode, float snapRatio, Frame fixedFrame, flan_CANCEL_ARG_CPP )
	{
	flan_FUNCTION_LOG;

	// Input validation
	if( source.isNull() ) return std::vector<Frame>();
	const Audio monoSource = source.convertToMono( canceller ); // Convert input to mono

	// Find expected waveform lengths via autocorrelation-based pitch detection
	const int acGranularity = 128;
	const int acWindowSize = 4096;
	std::vector<float> localWavelengths;
	Frame globalWavelength = 0;
	if( pitchMode != Wavetable::PitchMode::None )
		{
		localWavelengths = monoSource.getLocalWavelengths( 0, 0, -1, 2048, 128, nullptr, canceller );
		globalWavelength = monoSource.getAverageWavelength( localWavelengths, .5, 64 ); 
		if( globalWavelength == -1 ) 
			{
			std::cout << "Insignificant pitch data found in Wavetable construction.\n";
			return std::vector<Frame>();
			}
		}

	std::vector<Frame> waveformStarts;

	auto snapHandler = [snapMode, &monoSource, snapRatio, &canceller]( Frame frameToSnap, Frame snapSourceFrame, Frame maxSnap )
		{
		//const Frame maxSnap = snapRatio * std::abs( frameToSnap - snapSourceFrame );
		switch( snapMode )
			{
			case Wavetable::SnapMode::None: return frameToSnap;
			case Wavetable::SnapMode::Zero: return snapFrameToSample( monoSource, frameToSnap, 0, maxSnap, canceller );
			case Wavetable::SnapMode::Level: return snapFrameToSample( monoSource, frameToSnap, monoSource.getSample( 0, snapSourceFrame ), maxSnap, canceller );
			default: return frameToSnap;
			}
		};
	
	waveformStarts.push_back( snapHandler( 0, 0, snapRatio * globalWavelength ) );

	while( true ) // Eat until we run out of buffer
		{
		flan_CANCEL_POINT( std::vector<Frame>() );

		// Guess how many frames the next wavelength will be
		Frame expectedNumFrames = 0;
		if( pitchMode ==  Wavetable::PitchMode::Local )
			{
			const int localWavelengthIndex = std::floor( waveformStarts.back() / acGranularity );
			if( localWavelengthIndex >= localWavelengths.size() ) break;
			const Frame localWavelength_c = localWavelengths[localWavelengthIndex];
			if( localWavelength_c == 0 ) continue;

			expectedNumFrames =  localWavelength_c != -1? localWavelength_c : globalWavelength; 
			}
		else if( pitchMode == Wavetable::PitchMode::Global )
				expectedNumFrames = globalWavelength; 
		else if( pitchMode == Wavetable::PitchMode::None )
			expectedNumFrames = fixedFrame;

		if( waveformStarts.back() + expectedNumFrames >= monoSource.getNumFrames() ) break;

		// Snap expected waveform end if needed
		waveformStarts.push_back( snapHandler( waveformStarts.back() + expectedNumFrames, waveformStarts.back(), snapRatio * expectedNumFrames ) );
		}

    return waveformStarts;
	}

Wavetable::Wavetable( const Audio & source, SnapMode snapMode, PitchMode pitchMode, float snapRatio, Frame fixed, flan_CANCEL_ARG_CPP )
	: wavelength( 2048 )
    , table( resampleWaveforms( source, getWaveformStarts( source, snapMode, pitchMode, snapRatio, fixed, canceller ), wavelength, canceller ).setVolume( 1 ) )
	, numWaveforms( table.getNumFrames() / wavelength )
	{
	}

Wavetable::Wavetable( Func1x1 f, int numWaves, flan_CANCEL_ARG_CPP )
	: wavelength( 2048 )
    , table()
	, numWaveforms( numWaves )
	{
	const int oversample = 16;

	// Set up output
    Audio::Format format;
	format.numChannels = 1;
    format.numFrames = wavelength * numWaves;
	table = Audio( format );

	const int overLength = wavelength * oversample;
	std::vector<Sample> fBuffer( overLength * 3 );

	r8b::CDSPResampler rs( overLength, wavelength, fBuffer.size() );

	// For each waveform, resample to wavelength
	for( Waveform waveform = 0; waveform < numWaves; ++waveform )
		{
		flan_CANCEL_POINT();
		// Sample f into buffer
		for( int s = 0; s < overLength; ++s )
			{
			fBuffer[s] = f( waveform + float( s ) / overLength );
			fBuffer[s + 1 * wavelength * oversample] = fBuffer[s];
			fBuffer[s + 2 * wavelength * oversample] = fBuffer[s];
			}

		// Resample
		rs.oneshot( fBuffer.data(), fBuffer.size(), table.getSamplePointer( 0, wavelength * waveform ), wavelength );
		}
	}

Wavetable::Wavetable( const Audio & source, const std::vector<Frame> & waveformStarts, flan_CANCEL_ARG_CPP )
	: wavelength( 2048 )
    , table( resampleWaveforms( source, waveformStarts, wavelength, canceller ).setVolume( 1 ) )
	, numWaveforms( waveformStarts.size() - 1 )
	{
	}

Wavetable::Wavetable()
	: wavelength( 0 )
    , table( Audio() )
	, numWaveforms( 0 )
	{
	}

Audio Wavetable::synthesize( Time length, Func1x1 freq, Func1x1 index, Time granularityTime, flan_CANCEL_ARG_CPP ) const
    {
	// Ouput setup
	Audio::Format format = table.getFormat();
    format.numFrames = length * table.timeToFrame();
    Audio out( format );

	const Frame granularity = granularityTime * out.timeToFrame();

	// Resampler setup
	WDL_Resampler rs;
	rs.SetMode( true, 1, true, 512 );
	
	Frame outFramesGenerated = 0;
	Frame phaseFrame = 0;
	while( outFramesGenerated < out.getNumFrames() )
		{
		flan_CANCEL_POINT( Audio() );

		const double inFreq_c = double( table.getSampleRate() ) / wavelength;
		const double outFreq_c = freq( outFramesGenerated * out.frameToTime() );
		const float index_c = std::clamp( index( outFramesGenerated * out.frameToTime() ), 0.0f, float( numWaveforms - 1 ) );
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
			const Sample leftSample  = table.getSample( 0, ( phaseFrame + j ) % wavelength + wavelength * leftIndex  );
			const Sample rightsample = table.getSample( 0, ( phaseFrame + j ) % wavelength + wavelength * rightIndex );
			rsinbuf[j] = ( 1.0f - indexRemainder ) * leftSample + indexRemainder * rightsample;
			}
				
		// Resample directly into output buffer and update frame indices.
		outFramesGenerated += rs.ResampleOut( out.getSamplePointer( 0, outFramesGenerated ), wdlWanted, out.getNumFrames() - outFramesGenerated, 1 );
		phaseFrame = ( phaseFrame + wdlWanted ) % wavelength;
		}
	
	return out;
    }

void Wavetable::addFades( Frame fadeFrames, Waveform w )
	{
	auto singleTransform = [this, fadeFrames]( Waveform w )
		{
		Sample * waveform = getWaveform( w );
		for( Frame frame = 0; frame < fadeFrames - 1; ++frame )
			{
			const float fade = std::sin( pi / 2.0f * ( frame + 1 ) / fadeFrames );
			*(waveform + frame) *= fade;
			*(waveform + wavelength - 1 - frame) *= fade;
			}
		};

	if( w == -1 ) 
		for( Waveform i = 0; i < numWaveforms; ++i )
			singleTransform( i );
	else singleTransform( w );
	}

void Wavetable::removeJumps( Frame fadeFrames, Waveform w )
	{
	auto singleTransform = [this, fadeFrames]( Waveform w )
		{
		Sample * waveform = getWaveform( w );
		const Sample l = *waveform;
		const Sample r = *( waveform + wavelength - 1 );
		const Sample mid = ( l + r ) / 2.0f;

		auto fader = [mid]( Sample * s, float fade ){ *s = ( *s - mid ) * fade + mid; };

		for( Frame frame = 0; frame < fadeFrames - 1; ++frame )
			{
			const float fade = std::sin( pi / 2.0f * ( frame + 1 ) / fadeFrames );
			fader( waveform + frame, fade );
			fader( waveform + wavelength - 1 - frame, fade );
			}
		};

	if( w == -1 ) 
		for( Waveform i = 0; i < numWaveforms; ++i )
			singleTransform( i );
	else singleTransform( w );
	}

void Wavetable::removeDC( Waveform w )
	{
	auto singleTransform = [this]( Waveform w )
		{
		auto startIter = table.begin() + wavelength * w;
		const float dc = std::accumulate( startIter, startIter + wavelength, 0.0f ) / float( wavelength );
		std::for_each( startIter, startIter + wavelength, [dc]( Sample & s ){ s -= dc; } );
		};

	if( w == -1 ) 
		for( Waveform i = 0; i < numWaveforms; ++i )
			singleTransform( i );
	else singleTransform( w );
	}

Sample * Wavetable::getWaveform( int waveform )
	{
	return table.getSamplePointer( 0, waveform * wavelength );
	}
	