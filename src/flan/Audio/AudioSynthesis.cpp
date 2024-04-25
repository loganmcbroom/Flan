#include "flan/Audio/Audio.h"

#include <iostream>
#include <random>
#include <ctime>
#include <execution>
#include <ranges>
#include <algorithm>

#include "r8brain/CDSPResampler.h"
#include "WDL/resample.h"
#include "flan/WindowFunctions.h"
#include "flan/WindowFunctions.h"
#include "flan/FFTHelper.h"

#undef min
#undef max

namespace flan {

//=====================================================================================================================================
// Normal Synthesis
//=====================================================================================================================================

Audio Audio::synthesize_waveform( 
	const Function<Second, Amplitude> & wave, 
	Second length, 
	const Function<Second, Frequency> & freq, 
	FrameRate sample_rate, 
	int oversample 
	)
	{
	if( oversample < 1 || length <= 0 || sample_rate <= 0 )
		return Audio::create_null();

	// Set up output
	Audio::Format format;
	format.num_channels = 1;
	format.num_frames = length * sample_rate;
	format.sample_rate = sample_rate;
	Audio out( format );

	const Frame num_in_frames = out.get_num_frames() * oversample;
	const double in_sample_rate = sample_rate * oversample;

	auto frequency_sampled = freq.sample( 0, num_in_frames, 1.0f / in_sample_rate );

	// Get the phases that wave_samples needs to be evaluated at
	// Modular sum is associative, prefer exclusive_scan over partial_sum
	FunctionSample<float> phases = frequency_sampled.exclusive_scan( 0.0f, [&]( float a, float b ){ return std::fmod( a + b, in_sample_rate ); } );
	phases.for_each( [&]( float & phase ){ phase /= in_sample_rate; } ); // freq sum -> phase

	// Evaluate wave at phases
	std::vector<float> wave_samples( phases.size() );
	flan::for_each_i( wave_samples.size(), wave.get_execution_policy(), [&]( int i ){ wave_samples[i] = wave( phases[i] ); } );

	// Downsample to requested sample rate
	r8b::CDSPResampler resampler( in_sample_rate, out.get_sample_rate(), wave_samples.size() );
	resampler.oneshot<float, float>( &wave_samples[0], wave_samples.size(), &out.get_buffer()[0], out.get_buffer().size() );

	return out;
	}

Audio Audio::synthesize_white_noise(
	Second length,
	FrameRate sample_rate,
	int oversample 
	)
	{
	if( oversample < 1 || length <= 0 || sample_rate <= 0 )
		return Audio::create_null();
		
	std::random_device rd;
	std::mt19937 rng( rd() );
	std::uniform_real_distribution<> dis( -1.0f, 1.0f );

	Audio out = Audio::create_empty_with_length( length, 1, sample_rate * oversample );
	for( Frame frame = 0; frame < out.get_num_frames(); ++frame )
		out.get_sample( 0, frame ) = dis( rng );

	return out.resample( sample_rate );
	}

// This is based on the work of Phil Burk, as seen here: https://www.firstpr.com.au/dsp/pink-noise/phil_burk_19990905_patest_pink.c
Audio Audio::synthesize_pink_noise( 
	Second length,
	FrameRate sample_rate,
	int num_rows )
	{
	if( length <= 0 || sample_rate <= 0 || num_rows < 1 )
		return Audio::create_null();

	/* This converts a frame to a row in this pattern, with frame/row axes:
                     x  	
             x               x               
         x       x       x       x       
       x   x   x   x   x   x   x   x   
      x x x x x x x x x x x x x x x x 
	*/
	auto get_row = []( Frame index )
		{
		int num_zeros = 0;
		int n = index;
		while( (n & 1) == 0 )
			{
			n = n >> 1;
			num_zeros++;
			}
		return num_zeros;
		};

	std::random_device rd;
	std::mt19937 rng( rd() );
	std::uniform_real_distribution<> dis( -1.0f, 1.0f );

	Audio out = Audio::create_empty_with_length( length );

	float running_sum = 0;
	std::vector<float> rows( num_rows, 0 );
	for( Frame frame = 0; frame < out.get_num_frames(); ++frame )
		{
		const int index = frame % rows.size();

		if( index != 0 )
			{
			const int row = get_row( index );
			
			const float new_random = dis( rng );
			running_sum -= rows[row];
			running_sum += new_random;
			rows[row] = new_random;
			}
		
		const float sum = running_sum + dis( rng );

		out.get_sample( 0, frame ) = sum;
		}

	out.set_volume_in_place( 1 );

	return out;
	}

Audio Audio::synthesize_spectrum(
	Second length, 
	const Function<Second, Frequency> & freq,
	const Function<Harmonic, Frequency> & spread,
	const Function<Harmonic, Amplitude> & harmonic_scale,
	const Function<Frequency, Amplitude> & distribution,
	int fundamental_power,
	int spectrum_size_power,
	Channel num_channels,
	Second granularity_time
	)
	{
	if( length <= 0 ||
		fundamental_power <= 0 || 
		spectrum_size_power <= 0 ||
		fundamental_power > spectrum_size_power ||
		granularity_time <= 0 )
		return Audio::create_null(); 
	
	if( spectrum_size_power >= 32 )
		{
		std::cout << "I limited spectrum_size_power to that people wouldn't accidentally generate gigabytes of spectrum.";
		return Audio::create_null();
		}

	const Frequency fundamental = std::pow( 2, fundamental_power );
	const Frame wavelength = std::pow( 2, spectrum_size_power );

	Audio table = Audio::create_empty_with_frames( wavelength, num_channels );

	FFTHelper fft( wavelength, false, true, false );
	std::random_device rd;
	std::mt19937 rng( rd() );

	auto frequency_to_bin = [&]( Frequency f ) -> fBin { return f / table.get_sample_rate() * float( fft.complex_buffer_size() ); };
	auto bin_to_frequency = [&]( Bin b ) -> Frequency  { return b * table.get_sample_rate() / float( fft.complex_buffer_size() ); };

	const Harmonic num_harmonics = std::ceil( bin_to_frequency( fft.complex_buffer_size() ) / fundamental ) + 1;

	// Sample harmonic_scale
	std::vector<Amplitude> harmonic_scale_sampled( num_harmonics );
	flan::for_each_i( num_harmonics, harmonic_scale.get_execution_policy(), [&]( Harmonic h ){ harmonic_scale_sampled[h] = harmonic_scale(h+1); } );

	// Sample spread
	std::vector<Amplitude> spread_sampled( num_harmonics );
	flan::for_each_i( num_harmonics, spread.get_execution_policy(), [&]( Harmonic h ){ spread_sampled[h] = spread(h+1); } );

	flan::for_each_i( fft.complex_buffer_size(), distribution.get_execution_policy(), [&]( Bin bin )
		{
		const Frequency freq_c = bin_to_frequency( bin );
		const int harmonic = std::round( freq_c/fundamental );
		if( harmonic == 0 ) return;

		const Frequency harmonic_frequency = fundamental*harmonic;

		const Frequency epsilon = 0.001;
		const Frequency mean = harmonic_frequency;
		const Frequency sd = spread_sampled[ harmonic-1 ];
		const Frequency x = bin_to_frequency( bin );
		const Magnitude r = ( sd <= epsilon? x : distribution( (x - mean)/sd ) / sd ) * harmonic_scale_sampled[harmonic-1];

		const Radian theta = std::uniform_real_distribution( 0.0f, pi2 )( rng );

		fft.get_complex_buffer()[bin] = std::polar( r, theta );
		} );

	fft.c2r_execute();

	std::copy( fft.real_begin(), fft.real_begin() + table.get_num_frames(), table.get_buffer().begin() );

	//==============================================================================================================================================
	
	// Ouput setup
	Audio::Format format = table.get_format();
    format.num_frames = table.time_to_frame( length );
    Audio out( format );

	const Frame granularity = out.time_to_frame( granularity_time );

	flan::for_each_i( num_channels, freq.get_execution_policy(), [&]( Channel channel )
		{
		// Resampler setup
		WDL_Resampler rs;
		rs.SetMode( true, 1, true );

		const Frame channel_frame_jump = (float(channel)/num_channels) * wavelength;

		Frame out_frames_generated = 0;
		Frame phase_frame = 0;
		while( out_frames_generated < out.get_num_frames() )
			{
			const double in_freq_c = fundamental;
			const double out_freq_c = freq( out.frame_to_time( out_frames_generated ) );
			const double freq_scale = out_freq_c / in_freq_c;

			// Current resample setup
			rs.SetRates( out_freq_c, in_freq_c );
			WDL_ResampleSample * rsinbuf = nullptr;
			const Frame wdl_wanted = rs.ResamplePrepare( granularity, 1, &rsinbuf );

			// Copy waveform to wdl
			for( Frame j = 0; j < wdl_wanted; ++j )
				{
				const Sample sample = table.get_sample( 0, ( phase_frame + j + channel_frame_jump ) % wavelength );
				rsinbuf[j] = sample;
				}
					
			// Resample directly into output buffer and update frame indices.
			out_frames_generated += rs.ResampleOut( out.get_sample_pointer( channel, out_frames_generated ), 
				wdl_wanted, out.get_num_frames() - out_frames_generated, 1 );
			phase_frame = ( phase_frame + wdl_wanted ) % wavelength;
			}
		} );

	out.set_volume_in_place( 1 );
	
	return out;
	}

Audio Audio::synthesize_impulse(
	Frequency base_freq,
	Harmonic num_harmonics,
	float chroma,
	FrameRate sample_rate
	)
	{
	Frame num_frames = sample_rate / base_freq;
	if( num_frames % 2 == 0 ) ++num_frames;
	const Frame half_frames = (num_frames - 1) / 2; // Half, not including center

	float chroma_power = chroma;
	const float chroma_normalization = chroma == 1 ? 
		1.0f / num_harmonics : 
		( 1.0f - chroma ) / ( chroma - std::pow( chroma, num_harmonics + 1 ) );

	Audio impulse = Audio::create_empty_with_frames( num_frames, 1, sample_rate ); 
	for( Harmonic h = 1; h <= num_harmonics; ++h )
		{
		const Frequency harmonic_freq = h * base_freq;
		for( Frame frame = half_frames; frame < impulse.get_num_frames(); ++frame )
			{
			const Second time = impulse.frame_to_time( frame - half_frames );
			impulse.get_sample( 0, frame ) += chroma_power * chroma_normalization * std::cos( pi2 * harmonic_freq * time );
			}
		chroma_power *= chroma;
		}
		
	// Copy second half onto first half mirrored
	for( Frame frame = 0; frame < half_frames; ++frame )
		impulse.get_sample( 0, frame ) = impulse.get_sample( 0, num_frames - 1 - frame );

	return impulse;
	}

//=====================================================================================================================================
// Granular Controllers
//=====================================================================================================================================


std::vector<Second> integrate_event_rate( 
	Second length,
	const Function<Second, float> & events_per_second,
	const Function<Second, float> & scatter,
	FrameRate sample_rate
	)
	{
	const Frame length_frames = length * sample_rate;

	auto events_per_second_sampled = events_per_second.sample( 0, length_frames, 1.0f / sample_rate );
	events_per_second_sampled.for_each( []( float & eps ) { eps = std::max( 0.0f, eps ); } );

	auto scatter_sampled = scatter.sample( 0, length_frames, 1.0f / sample_rate );
	scatter_sampled.for_each( []( float & scatter ) { scatter = std::max( 0.0f, scatter ); } );

	// Generate unscattered events by integrating eventsPerFrame
	// When the integral passes an integer, an event occurs
	std::vector<Frame> event_frames;
	float event_accumulator = 1.0f;
	for( Frame frame = 0; frame < length_frames; ++frame )
		{
		const Second t = fFrame( frame ) / sample_rate;
		event_accumulator += events_per_second_sampled[frame] / sample_rate;
		if( event_accumulator >= 1.0f )
			{
			event_frames.push_back( frame );
			event_accumulator -= std::floor( event_accumulator );
			}
		}
	
	// Scatter event frames, removing those landing outside length
	// This can cause the events around the boundaries to be less dense than expected
	std::default_random_engine rng( std::time( nullptr ) );
	std::for_each( event_frames.begin(), event_frames.end(), [&]( Frame & event_frame )
		{
		// Scatter is standard deviation in seconds.
		const Second t = fFrame( event_frame ) / sample_rate;
		const float scatter_t = scatter_sampled[event_frame];
		const float events_per_second_t = events_per_second_sampled[event_frame];
		if( scatter_t == 0 ) return; // If there is no scatter, there is nothing to do
		if( events_per_second_t == 0 ) return; // If there is no EpS, there is nothing to do

		const Second std_in_seconds = scatter_t / events_per_second_t;
		const fFrame std_in_frames = std_in_seconds * sample_rate;
		const Frame scattered_frame = std::normal_distribution<fFrame>( event_frame, std_in_frames )( rng );
		if( 0 <= scattered_frame && scattered_frame < length_frames )
			event_frame = scattered_frame;
		else event_frame = std::numeric_limits<Frame>::max(); // Using max to indicate "don't use", it will be sorted to the end shortly
		} );

	// Get our ducks back in a row
	std::sort( FLAN_PAR_UNSEQ event_frames.begin(), event_frames.end() ); 

	// Remove anything set to max at the end of the events; these were scattered outside the bounds
	auto oob_event_iter = std::lower_bound( event_frames.begin(), event_frames.end(), std::numeric_limits<Frame>::max() );
	if( oob_event_iter != event_frames.end() )
		event_frames.erase( oob_event_iter, event_frames.end() );

	// Convert event frames back to times
	std::vector<Second> event_times( event_frames.size() );
	std::transform( FLAN_PAR_UNSEQ event_frames.begin(), event_frames.end(), event_times.begin(), 
		[&]( Frame frame ){ return frame / sample_rate; } );

	return event_times;
	}

Audio Audio::synthesize_grains( 
	Second length,
	const Function<Second, float> & grains_per_second,
	const Function<Second, float> & time_scatter,
	const Function<Second, Audio> & grain_source,
	FrameRate sample_rate
	)
	{
	// Input validation
	if( length <= 0 ) return Audio::create_null();
	
	const std::vector<Second> event_times = integrate_event_rate( length, grains_per_second, time_scatter, sample_rate );

	// Generate trainlets
	std::vector<Audio> grains( event_times.size() );
	flan::for_each_i( event_times.size(), grain_source.get_execution_policy(), [&]( Index i )
		{
		const Second t = event_times[i];
		grains[i] = grain_source( t );
		} );
	
	return Audio::mix( grains, event_times );
	}

// Synth_grains optimized for single input
Audio synthesize_grains_repeat(
	const Audio & me,
	Second length,
	const Function<Second, float> & grains_per_second,
	const Function<Second, float> & time_scatter,
	const Function<Second, Amplitude> & gain
	)
	{
	if( me.is_null() ) return Audio::create_null();

	// Input validation
	if( length <= 0 ) return Audio::create_null();
	
	const std::vector<Second> event_times = integrate_event_rate( length, grains_per_second, time_scatter, me.get_sample_rate() );

	std::vector<Amplitude> gains( event_times.size() );
	for( int i = 0; i < event_times.size(); ++i )
		gains[i] = gain( event_times[i] );

	return Audio::mix( std::vector<const Audio *>( event_times.size(), &me ), event_times, gains );
	}

Audio Audio::texture( 
	Second length, 
	const Function<Second, float> & grains_per_second, 
	const Function<Second, float> & time_scatter, 
	const AudioMod & mod,
	bool mod_feedback
	) const
	{
	if( is_null() ) return Audio::create_null();

	if( mod.is_null() )
		return synthesize_grains_repeat( *this, length, grains_per_second, time_scatter, 1.0f );

	const std::vector<Second> event_times = integrate_event_rate( length, grains_per_second, time_scatter, get_sample_rate() );

	std::vector<Audio> modded_audio;
	modded_audio.reserve( event_times.size() );
	if( mod_feedback )
		{
		Audio this_copy = copy();
		mod( this_copy, 0 );
		modded_audio.emplace_back( std::move( this_copy ) );
		// Sequential execution policy is required by feedback
		std::transform( event_times.begin() + 1, event_times.end(), std::back_inserter( modded_audio ), [&]( Second event_time )
			{
			Audio input = modded_audio.back().copy();
			mod( input, event_time );
			return input;
			} );
		}
	else // Convert this using mod at each event time according to the mod execution policy
		{ // This could be handled by synthesize_grains, but it would make things a little harder in composition functions
		modded_audio.resize( event_times.size() );
		flan::for_each_i( event_times.size(), mod.get_execution_policy(), [&]( Index i )
			{
			const Second event_time = event_times[i];
			Audio this_copy = copy();
			mod( this_copy, event_time );
			modded_audio[i] = std::move( this_copy );
			} );
		// runtime_execution_policy_handler( mod.get_execution_policy(), [&]( auto policy ) {
		// 	std::transform( FLAN_POLICY event_times.begin(), event_times.end(), modded_audio.begin(), [&]( const Second event_time ) 
		// 		{
		// 		Audio this_copy = copy();
		// 		mod( this_copy, event_time );
		// 		return std::move( this_copy );
		// 		} ); 
		// 	} );
		}
	return Audio::mix( modded_audio, event_times );
	}

//=====================================================================================================================================
// Granular Compositions
//=====================================================================================================================================

Audio Audio::synthesize_trainlets( 
	Second length,
	const Function<Second, float> & grains_per_second,
	const Function<Second, float> & time_scatter,
	const Function<Second, vec2> & position,
	const Function<Second, Amplitude> & trainlet_gain_envelope,
	const Function<Second, Frequency> & freq,
	const Function<Second, Second> & trainlet_length,
	const Function<Second, Harmonic> & num_harmonics,
	const Function<Second, float> & chroma,
	const Function<Second, Frequency> & impulse_harmonic_frequency,
	FrameRate sample_rate
	)
	{
	const ExecutionPolicy ex_pol = lowest_execution( position, trainlet_gain_envelope, freq, trainlet_length, num_harmonics, chroma, impulse_harmonic_frequency );

	return synthesize_grains( length, grains_per_second, time_scatter, Function<Second, Audio>( [&]( Second t )
		{ 
		const Audio impulse = synthesize_impulse( impulse_harmonic_frequency( t ), num_harmonics( t ), chroma( t ), sample_rate );
		return synthesize_grains_repeat(
			impulse,
			trainlet_length( t ),
			freq,
			0,
			trainlet_gain_envelope
			).stereo_spatialize( position( t ) );
		}, ex_pol ) );
	}

Audio Audio::granulate( 
	Second length, 
	const Function<Second, float> & grains_per_second, 
	const Function<Second, float> & time_scatter, 
	const Function<Second, Second> & time_selection, 
	const Function<Second, Second> & grain_length,
	const Function<Second, Second> & fade_time,
	const AudioMod & mod
	) const
	{
	if( is_null() ) return Audio::create_null();

	const auto selection_sampled 	= time_selection.sample( 0, time_to_frame( length ), frame_to_time( 1 ) );
	const auto grain_length_sampled = grain_length	.sample( 0, time_to_frame( length ), frame_to_time( 1 ) );
	const auto fade_time_sampled 	= fade_time		.sample( 0, time_to_frame( length ), frame_to_time( 1 ) );

	auto grain_generator = [&]( const Audio & in, float t ) -> Audio
		{
		const Second selection_c = selection_sampled[ time_to_frame( t ) ];
		const Second grain_length_c = grain_length_sampled[ time_to_frame( t ) ];
		const Second fade_time_c = fade_time_sampled[ time_to_frame( t ) ];
		Audio grain = in.cut( selection_c, selection_c + grain_length_c, fade_time_c, fade_time_c );
		return std::move( grain );
		};

	return synthesize_grains( 
		length, 
		grains_per_second, 
		time_scatter,
		Function<Second, Audio>( [&]( Second t )
			{ 
			Audio grain = grain_generator( *this, t );
			if( ! mod.is_null() )
				mod( grain, t );
			return std::move( grain );
			}, mod.get_execution_policy() ),
		get_sample_rate() );
	}

Audio Audio::psola(
	Second length, 
	const Function<Second, Second> & time_selection,
	const AudioMod & mod
	) const
	{
	auto freq = get_frequency_envelope();

	// A bit janky to sample here and wrap the samples in a lambda, but it's an easy way to protect from double calling
	auto time_selection_sampled = time_selection.sample( 0, std::ceil( time_to_frame( length ) ), frame_to_time( 1 ) );
	AudioMod composition_mod = [&]( Audio & a, Second t )
		{
		if( ! mod.is_null() )
			mod( a, t );

		a.modify_volume_in_place( [&a]( Second t ){ return Windows::hann( t / a.get_length() ); } );
		};

	return granulate( 
		length,
		[&]( Second t ){ return freq( time_selection_sampled[time_to_frame(t)] ); },
		0,
		[&]( Second t ){ return time_selection_sampled[time_to_frame(t)]; },
		[&]( Second t ){ return 2.0f / freq( time_selection_sampled[time_to_frame(t)] ); },
		0.05,
		composition_mod 
		);
	}

}
