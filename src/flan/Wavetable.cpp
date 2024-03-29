#include "flan/Wavetable.h"

#include <algorithm>
#include <numeric>
#include <iostream>
//#include <ranges>

#include "WDL/resample.h"
#include "r8brain/CDSPResampler.h"
#include "flan/PV/PV.h"
#include "flan/Graph.h"
#include "flan/FFTHelper.h"

#undef min
#undef max

using namespace flan;

static Frame snap_frame_to_sample( const Sample * data, int n, Frame frame_to_snap, Sample snapHeight, Frame search_distance )
	{	
	const Frame left_frame_limit = std::max( frame_to_snap - search_distance, 0 );
	const Frame right_frame_limit = std::min( frame_to_snap + search_distance, n - 1 );

	// Bidirectionally search from start_frame until a crossing is found
	const bool is_above = data[frame_to_snap] > snapHeight;
	for( Frame search_offset = 0; search_offset <= search_distance; ++search_offset )
		{
		// Check left
		const Frame left_search_frame = frame_to_snap - search_offset;
		if( left_search_frame >= left_frame_limit && data[left_search_frame] > snapHeight != is_above )
			return left_search_frame + 1;
			
		// Check right
		const Frame right_search_frame = frame_to_snap + search_offset;
		if( right_search_frame < right_frame_limit && data[right_search_frame] > snapHeight != is_above )
			return right_search_frame;
		}
	
	// If control reaches here, cross search failed. Use frame with output nearest to crossing.
	auto norm = [&]( Frame f )
		{ 
		const float m = 1.0f;
		const float r = 1.0f + m * float( std::abs( f - frame_to_snap ) ) / search_distance; // Slightly prefer frames nearer to initial
		return std::abs( data[f] - snapHeight ) * r;
		};

	Frame current_best_frame = frame_to_snap;
	Sample current_best_sample_distance = norm( frame_to_snap );
	for( Frame search_frame = left_frame_limit; search_frame <= right_frame_limit; ++search_frame )
		{
		const Sample current_sample_distance = norm( search_frame );
		if( current_sample_distance < current_best_sample_distance )
			{
			current_best_frame = search_frame;
			current_best_sample_distance = current_sample_distance;
			}
		}

	return current_best_frame;
	}

static Frame snap_frame_to_sample( const Audio & me, Channel channel, Frame frame_to_snap, Sample snapHeight, Frame search_distance )
	{
	return snap_frame_to_sample( me.get_sample_pointer( channel, 0 ), me.get_num_frames(), frame_to_snap, snapHeight, search_distance );
	}

static Audio resample_waveforms( const Audio & source, const std::vector<std::vector<Frame>> & waveform_starts, Frame output_wavelength )
	{
	// Input validation
	if( source.is_null() ) return Audio::create_null();
	if( waveform_starts.empty() ) return Audio::create_null();
	for( const auto & w : waveform_starts ) 
		if( w.empty() ) return Audio::create_null();

	// Set up output
    Audio::Format format = source.get_format();
	const size_t max_num_starts = std::max_element( waveform_starts.begin(), waveform_starts.end(), []( const std::vector<Frame> & a, const std::vector<Frame> & b )
		{ return a.size() < b.size(); } )->size();
    format.num_frames = output_wavelength * max_num_starts;
	Audio out( format );
	out.clear_buffer();

	FFTHelper ifft( output_wavelength, false, true, true );

	for( Channel channel = 0; channel < source.get_num_channels(); ++channel )
		{	
		// For each waveform, resample to wavelength
		for( int waveform = 0; waveform < waveform_starts[channel].size() - 1; ++waveform )
			{
			const Frame start_frame = waveform_starts[channel][waveform];
			const Frame end_frame =  waveform_starts[channel][waveform + 1];
			const Frame num_input_frames = end_frame - start_frame;
			const float expansion = float( output_wavelength ) / num_input_frames;

			FFTHelper fft( num_input_frames, true, false, false );
			std::copy( source.channel_begin( channel ) + start_frame, source.channel_begin( channel ) + end_frame, fft.real_begin() );
			fft.r2c_execute();

			std::fill( ifft.complex_begin(), ifft.complex_end(), 0 );
			std::copy( fft.complex_begin(), fft.complex_end(), ifft.complex_begin() );
			ifft.c2r_execute();

			// We would love to just copy directly to out, but the resampling process will usually shift the zero-crossing, so we need to realign
			// the output so it sits at the start of the waveform.
			Frame zero_crossing_frame = 0;
			const Frame search_distance = output_wavelength * 0.1f;
			const bool is_above = ifft.get_real_buffer()[0] > 0;
			for( Frame search_offset = 1; search_offset <= search_distance; ++search_offset )
				{
				if( ifft.get_real_buffer()[output_wavelength - search_offset] > 0 != is_above )
					{
					zero_crossing_frame = output_wavelength - search_offset;
					break;
					}
				if( ifft.get_real_buffer()[search_offset] > 0 != is_above )
					{
					zero_crossing_frame = search_offset;
					break;
					}
				}

			for( Frame f = 0; f < output_wavelength; ++f )
				{
				const Sample unscaled_output = ifft.get_real_buffer()[( zero_crossing_frame + f ) % output_wavelength];
				out.get_sample( channel, waveform * output_wavelength + f ) = unscaled_output / std::sqrt( num_input_frames * num_input_frames );
				}
			};
		}

	//out.set_volume_in_place( 1 );
    return out;
    }

static std::vector<std::vector<Frame>> get_waveform_starts( 
	const Audio & source, 
	Wavetable::SnapMode snap_mode, 
	Wavetable::PitchMode pitch_mode, 
	Frame wavelength,
	float snap_ratio, 
	Frame fixed_frame, 
	flan_CANCEL_ARG_CPP )
	{
	// Input validation
	if( source.is_null() ) return std::vector<std::vector<Frame>>();
	if( fixed_frame < 1 )  return std::vector<std::vector<Frame>>();
	if( snap_ratio <= 0 || snap_ratio >= .95 ) return std::vector<std::vector<Frame>>();

	const Audio lp_source = source.filter_1pole_lowpass( 4000, 2 );

	std::vector<std::vector<Frame>> waveform_starts( source.get_num_channels() );

	// Find expected waveform lengths via autocorrelation-based pitch detection
	const Frame ac_granularity = 128;
	const Frame ac_windowSize = wavelength * 2;
	for( Channel channel = 0; channel < source.get_num_channels(); ++channel )
		{
		std::vector<fFrame> local_wavelengths;
		Frame global_wavelength = 0;
		if( pitch_mode != Wavetable::PitchMode::None )
			{
			local_wavelengths = lp_source.get_local_wavelengths( channel, 0, -1, wavelength, ac_granularity, 1, 32, canceller );
			global_wavelength = lp_source.get_average_wavelength( local_wavelengths, .2, 64 ); 
			if( pitch_mode == Wavetable::PitchMode::Global && global_wavelength == -1 ) 
				pitch_mode = Wavetable::PitchMode::None;
			}
		
		auto snap_handler = [&]( Frame frame_to_snap, Frame snapSourceFrame, Frame max_snap )
			{
			switch( snap_mode )
				{
				case Wavetable::SnapMode::None: 
					return frame_to_snap;
				case Wavetable::SnapMode::Zero: 
					return snap_frame_to_sample( source, channel, frame_to_snap, 0, max_snap );
				case Wavetable::SnapMode::Level: 
					return snap_frame_to_sample( source, channel, frame_to_snap, source.get_sample( channel, snapSourceFrame ), max_snap );
				default: return frame_to_snap;
				}
			};

		waveform_starts[channel].push_back( snap_handler( 0, 0, snap_ratio * global_wavelength ) );
		while( true ) // Eat until we run out of buffer
			{
			if( canceller ) return std::vector<std::vector<Frame>>();

			// Guess how many frames the next wavelength will be
			Frame expected_num_frames = 0;
			if( pitch_mode ==  Wavetable::PitchMode::Local )
				{
				// Where in the local wavelegnths vector we currently are
				const int local_wavelength_index = std::floor( waveform_starts[channel].back() / ac_granularity );
				if( local_wavelength_index >= local_wavelengths.size() ) break;
				const Frame local_wavelength_c = local_wavelengths[local_wavelength_index];

				if( local_wavelength_c > 0 )
					expected_num_frames = local_wavelength_c;
				else if( global_wavelength > 0 )
					expected_num_frames = global_wavelength;
				else 
					expected_num_frames = fixed_frame;
				}
			else if( pitch_mode == Wavetable::PitchMode::Global )
					expected_num_frames = global_wavelength; 
			else if( pitch_mode == Wavetable::PitchMode::None )
				expected_num_frames = fixed_frame;

			if( waveform_starts[channel].back() + expected_num_frames >= source.get_num_frames() ) break;

			// Snap expected waveform end if needed
			waveform_starts[channel].push_back( snap_handler( 
				waveform_starts[channel].back() + expected_num_frames, 
				waveform_starts[channel].back(), 
				snap_ratio * expected_num_frames ) );
			}
		}

    return waveform_starts;
	}

Wavetable::Wavetable( 
	const Audio & source, 
	SnapMode snap_mode, 
	PitchMode pitch_mode, 
	Frame _wavelength,
	float snap_ratio, 
	Frame fixed, 
	flan_CANCEL_ARG_CPP )
	: wavelength( _wavelength )
	, num_source_frames( source.get_num_frames() )
	, waveform_starts( get_waveform_starts( source, snap_mode, pitch_mode, wavelength, snap_ratio, fixed, canceller ) )
    , table( resample_waveforms( source, waveform_starts, wavelength ) )
	{
	}

Wavetable::Wavetable( const Function<Second, Amplitude> & f, int num_waves, Frame _wavelength, flan_CANCEL_ARG_CPP )
	: wavelength( _wavelength )
	, num_source_frames( num_waves )
	, waveform_starts( 1, std::vector<Frame>( num_waves ) )
    , table( Audio::create_empty_with_frames( wavelength * num_waves, 1, 48000 ) )
	{
	for( int i = 0; i < num_waves; ++i )
		waveform_starts[0][i] = i; 

	// For each waveform, resample to wavelength
	for( Waveform waveform = 0; waveform < num_waves; ++waveform )
		for( Frame frame = 0; frame < wavelength; ++frame )
			table.get_sample( 0, frame ) = f( waveform + float( frame ) / wavelength );
	}

// Wavetable::Wavetable( const Audio & source, const std::vector<std::vector<Frame>> & waveform_starts, flan_CANCEL_ARG_CPP )
// 	: wavelength( _wavelength )
//     , table()
// 	, table_info()
// 	{
// 	resample_waveforms( source, waveform_starts, wavelength, canceller );
// 	}

Wavetable::Wavetable()
	: wavelength( 0 )
	, num_source_frames( 0 )
	, waveform_starts()
    , table( Audio::create_null() )
	{
	}

Audio Wavetable::synthesize( 
	Second length, 
	const Function<Second, Frequency> & freq, 
	const Function<Second, float> & ratio, 
	bool smooth,
	Second granularity_time, 
	flan_CANCEL_ARG_CPP 
	) const
    {
	if( is_null() ) return Audio::create_null();
	
	// Ouput setup
	Audio::Format format = table.get_format();
    format.num_frames = table.time_to_frame( length );
    Audio out( format );

	const Frame granularity = out.time_to_frame( granularity_time );

	const ExecutionPolicy lowest_policy = lowest_execution( freq, ratio );

	//for( Channel channel = 0; channel < table.get_num_channels(); ++channel )
	flan::for_each_i( table.get_num_channels(), lowest_policy, [&]( Channel channel )
		{
		// Resampler setup
		WDL_Resampler rs;
		rs.SetMode( true, 1, true );

		Frame out_frames_generated = 0;
		Frame phase_frame = 0;
		while( out_frames_generated < out.get_num_frames() )
			{
			if( canceller ) return;

			const double in_freq_c = double( table.get_sample_rate() ) / wavelength;
			const double out_freq_c = freq( out.frame_to_time( out_frames_generated ) );
			const float table_index_c = ratio_to_table_index( ratio( out.frame_to_time( out_frames_generated ) ), channel );
			const Frame left_index = std::floor( table_index_c );
			const Frame right_index = std::ceil( table_index_c );
			const float index_remainder = table_index_c - left_index;

			// Current resample setup
			rs.SetRates( out_freq_c, in_freq_c );
			WDL_ResampleSample * rsinbuf = nullptr;
			const Frame wdl_wanted = rs.ResamplePrepare( granularity, 1, &rsinbuf );

			// Copy waveform to wdl
			for( Frame j = 0; j < wdl_wanted; ++j )
				{
				if( smooth )
					{
					const Sample left_sample  = table.get_sample( channel, ( phase_frame + j ) % wavelength + wavelength * left_index  );
					const Sample right_sample = table.get_sample( channel, ( phase_frame + j ) % wavelength + wavelength * right_index );
					rsinbuf[j] = ( 1.0f - index_remainder ) * left_sample + index_remainder * right_sample;
					}
				else
					{
					const Sample left_sample  = table.get_sample( channel, ( phase_frame + j ) % wavelength + wavelength * left_index  );
					rsinbuf[j] = left_sample;
					}
				}
					
			// Resample directly into output buffer and update frame indices.
			out_frames_generated += rs.ResampleOut( out.get_sample_pointer( channel, out_frames_generated ), wdl_wanted, out.get_num_frames() - out_frames_generated, 1 );
			phase_frame = ( phase_frame + wdl_wanted ) % wavelength;
			}
		} );
	
	return out;
    }

Graph Wavetable::graph_waveform_range( Channel channel, int start, int num ) const
	{
	if( is_null() ) return Graph();

	Graph g( -1, -1 );
	g.fill_image( Color::from_hsv( 0, 0, .04 ) );
	g.add_full_split_view_y( Rect( 0, 0, 1, 1 ), num );
		
	std::vector<const float *> channel_starts;
	for( int i = start; i < start + num; ++i )
		if( 0 <= i && i < get_num_waveforms( channel ) )
			channel_starts.push_back( get_waveform( i, channel ) );
	g.draw_waveforms( channel_starts, wavelength );

	return g;
	}

void Wavetable::save_waveform_range_to_bmp( const std::string & filename, Channel channel, int start, int end ) const
	{
	const Graph g = graph_waveform_range( channel, start, end );
	g.save_image( filename );
	}

int Wavetable::get_num_waveforms( Channel c ) const
	{
	return waveform_starts[c].size();
	}

void Wavetable::add_fades_in_place( Frame fade_frames )
	{
	if( is_null() ) return;

	flan::for_each_i( waveform_starts.size(), ExecutionPolicy::Parallel_Unsequenced, [&]( Channel channel )
		{
		flan::for_each_i( waveform_starts[channel].size(), ExecutionPolicy::Parallel_Unsequenced, [&]( int waveform_index )
			{
			for( Channel channel = 0; channel < table.get_num_channels(); ++channel )
				{
				Sample * waveform = get_waveform( waveform_index, channel );
				for( Frame frame = 0; frame < fade_frames - 1; ++frame )
					{
					const float fade = std::sin( pi / 2.0f * ( frame + 1 ) / fade_frames );
					*(waveform + frame) *= fade;
					*(waveform + wavelength - 1 - frame) *= fade;
					}
				}
			} );
		} );
	}

void Wavetable::remove_jumps_in_place( Frame fade_frames )
	{
	if( is_null() ) return;

	flan::for_each_i( waveform_starts.size(), ExecutionPolicy::Parallel_Unsequenced, [&]( Channel channel )
		{
		flan::for_each_i( waveform_starts[channel].size(), ExecutionPolicy::Parallel_Unsequenced, [&]( int waveform_index )
			{
			for( Channel channel = 0; channel < table.get_num_channels(); ++channel )
				{
				Sample * waveform = get_waveform( waveform_index, channel );
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
				}
			} );
		} );
	}

void Wavetable::remove_dc_in_place()
	{
	if( is_null() ) return;

	flan::for_each_i( waveform_starts.size(), ExecutionPolicy::Parallel_Unsequenced, [&]( Channel channel )
		{
		flan::for_each_i( waveform_starts[channel].size(), ExecutionPolicy::Parallel_Unsequenced, [&]( int waveform_index )
			{
			Sample * start = get_waveform( waveform_index, channel );
			const float dc = std::accumulate( start, start + wavelength, 0.0f ) / float( wavelength );
			std::for_each( FLAN_PAR_UNSEQ start, start + wavelength, [&]( Sample & s ){ s -= dc; } );
			} );
		} );
	}

void Wavetable::normalize_in_place()
	{
	if( is_null() ) return;

	flan::for_each_i( waveform_starts.size(), ExecutionPolicy::Parallel_Unsequenced, [&]( Channel channel )
		{
		flan::for_each_i( waveform_starts[channel].size(), ExecutionPolicy::Parallel_Unsequenced, [&]( int waveform_index )
			{
			Sample * start = get_waveform( waveform_index, channel );

			Magnitude max_amp = 0;
			for( Frame frame = 0; frame < wavelength; ++frame )
				{
				const Sample & s = *(start + frame);
				if( std::abs( s ) > max_amp )
					max_amp = std::abs( s );
				}

			if( max_amp < 0.001 ) return;
			std::for_each( FLAN_PAR_UNSEQ start, start + wavelength, [&]( Sample & s ){ s /= max_amp; } );
			} );
		} );
	}

Sample * Wavetable::get_waveform( int waveform_index, Channel channel )
	{
	return table.get_sample_pointer( channel, waveform_index * wavelength );
	}

const Sample * Wavetable::get_waveform( int waveform_index, Channel channel ) const
	{
	return table.get_sample_pointer( channel, waveform_index * wavelength );
	}

float Wavetable::ratio_to_table_index( float r, Channel channel ) const
	{
	/*
	Waveforms are not usually sourced from fixed-size sections of Audio, rather the domains they occupied in time in the source vary.
	We would like to still maintain temporal behaviour matching the source as we scroll through the table.
	This conversion handles that.
	*/

	const Frame source_frame = r * num_source_frames;

	if( source_frame <= 0 ) return 0;
	if( source_frame > num_source_frames ) return waveform_starts[channel].size() - 1;

	const auto right_iter = std::upper_bound( waveform_starts[channel].begin(), waveform_starts[channel].end(), source_frame );
	const int right_index = std::distance( waveform_starts[channel].begin(), right_iter );

	if( right_index == 0 ) return 0;
	if( right_index == waveform_starts[channel].size() ) return waveform_starts[channel].size() - 1;

	const Frame left_source_frame = *(right_iter - 1);
	const Frame right_source_frame = *right_iter;

	const float index = right_index - 1 + float( source_frame - left_source_frame ) / ( right_source_frame - left_source_frame );

	return std::clamp( index, 0.0f, float( waveform_starts[channel].size() - 1 ) );
	}

bool Wavetable::is_null() const
	{
	if( wavelength <= 0 ) return true;

	if( waveform_starts.empty() ) return true;
	for( Channel channel = 0; channel < waveform_starts.size(); ++channel )
		if( waveform_starts[channel].empty() ) return true;

	if( num_source_frames <= 0 ) return true;

	if( table.is_null() ) return true;
    
	return false;
	}
	