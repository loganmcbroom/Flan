#include "flan/Audio/Audio.h"

#include <iostream>
#include <random>
#include <ctime>
#include <execution>
#include <ranges>
#include <algorithm>

#include "r8brain/CDSPResampler.h"
#include "flan/WindowFunctions.h"
#include "flan/WindowFunctions.h"

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
	size_t samplerate, 
	size_t oversample 
	)
	{
	if( length <= 0 ) return Audio::create_null();

	// Set up output
	Audio::Format format;
	format.num_channels = 1;
	format.num_frames = length * samplerate;
	format.sample_rate = samplerate;
	Audio out( format );

	const Frame num_in_frames = out.get_num_frames() * oversample;
	const double in_sample_rate = samplerate * oversample;

	auto frequency_sampled = freq.sample( 0, num_in_frames, 1.0f / in_sample_rate );

	// Get the phases that wave_samples needs to be evaluated at
	std::vector<float> phases( frequency_sampled.size() );
	std::exclusive_scan( FLAN_PAR_UNSEQ frequency_sampled.begin(), frequency_sampled.end(), phases.begin(), 0.0f, [&]( float a, float b ) -> float { 
		return std::fmod( a + b, in_sample_rate ); } ); // Modular sum is associative, prefer exclusive_scan over partial_sum
	std::for_each( FLAN_PAR_UNSEQ phases.begin(), phases.end(), [&]( float & phase ){ phase /= in_sample_rate; } ); // freq sum -> phase

	// Evaluate wave at phases using wave execution policy
	std::vector<float> wave_samples( phases.size() );
	runtime_execution_policy_handler( wave.get_execution_policy(), [&]( auto policy ){
		std::transform( FLAN_POLICY phases.begin(), phases.end(), wave_samples.begin(), [&]( const float phase ){ return wave( phase ); } );
		} );

	// Downsample to requested sample rate
	r8b::CDSPResampler resampler( in_sample_rate, out.get_sample_rate(), wave_samples.size() );
	resampler.oneshot<float, float>( &wave_samples[0], wave_samples.size(), &out.get_buffer()[0], out.get_buffer().size() );

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
	std::for_each( events_per_second_sampled.begin(), events_per_second_sampled.end(), []( float & eps )
		{ 
		eps = std::max( 0.0f, eps );
		} );

	auto scatter_sampled = scatter.sample( 0, length_frames, 1.0f / sample_rate );
	std::for_each( scatter_sampled.begin(), scatter_sampled.end(), []( float & scatter )
		{ 
		scatter = std::max( 0.0f, scatter );
		} );

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
		const Frame scattered_frame = std::normal_distribution<Second>( event_frame, std_in_frames )( rng );
		if( 0 <= scattered_frame && scattered_frame < length_frames )
			event_frame = scattered_frame;
		else event_frame = std::numeric_limits<Frame>::max(); // Using max to indicate "don't use", it will be sorted to the end shortly
		} );

	// Get our ducks back in a row
	std::sort( FLAN_PAR_UNSEQ event_frames.begin(), event_frames.end() ); 

	// Remove anything set to max at the end of the events; these were scattered outside the bounds
	auto oob_event_iter = std::lower_bound( event_frames.begin(), event_frames.end(), std::numeric_limits<Frame>::max() );
	if( oob_event_iter != event_frames.end() )
		event_frames.erase( oob_event_iter );

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
	std::cout << __FUNCTION__ << std::endl;

	// Input validation
	if( length <= 0 ) return Audio::create_null();
	
	const std::vector<Second> event_times = integrate_event_rate( length, grains_per_second, time_scatter, sample_rate );

	// Generate trainlets
	std::vector<Audio> grains( event_times.size() );
	runtime_execution_policy_handler( grain_source.get_execution_policy(), [&]( auto policy )
		{
		std::for_each( FLAN_POLICY iota_iter( 0 ), iota_iter( event_times.size() ), [&]( Index i )
			{
			const Second t = event_times[i];
			grains[i] = grain_source( t );
			} );
		} );
	
	return Audio::mix( grains, event_times );
	}

Audio Audio::synthesize_grains_repeat(
	Second length,
	const Function<Second, float> & grains_per_second,
	const Function<Second, float> & time_scatter,
	const Function<Second, Amplitude> & gain,
	FrameRate sample_rate
	) const
	{
	if( is_null() ) return Audio::create_null();

	// Input validation
	if( length <= 0 ) return Audio::create_null();
	
	const std::vector<Second> event_times = integrate_event_rate( length, grains_per_second, time_scatter, sample_rate );

	std::vector<Amplitude> gains( event_times.size() );
	for( int i = 0; i < event_times.size(); ++i )
		gains[i] = gain( event_times[i] );

	return Audio::mix( std::vector<const Audio *>( event_times.size(), this ), event_times, gains );
	}

Audio Audio::synthesize_grains_with_feedback_mod( 
	Second length, 
	const Function<Second, float> & grains_per_second, 
	const Function<Second, float> & time_scatter, 
	const AudioMod & mod,
	bool mod_feedback, 
	FrameRate sample_rate
	) const
	{
	if( is_null() ) return Audio::create_null();

	if( mod.is_null() )
		return synthesize_grains_repeat( length, grains_per_second, time_scatter, 1.0f, sample_rate );

	const std::vector<Second> event_times = integrate_event_rate( length, grains_per_second, time_scatter, sample_rate );

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
		runtime_execution_policy_handler( mod.get_execution_policy(), [&]( auto policy ) {
			std::transform( FLAN_POLICY event_times.begin(), event_times.end(), modded_audio.begin(), [&]( const Second event_time ) 
				{
				Audio this_copy = copy();
				mod( this_copy, event_time );
				return std::move( this_copy );
				} ); 
			} );
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
	std::cout << __FUNCTION__ << std::endl;

	const ExecutionPolicy ex_pol = lowest_execution( 
			position.get_execution_policy(), 
			trainlet_gain_envelope.get_execution_policy(),
			freq.get_execution_policy(),
			trainlet_length.get_execution_policy(),
			num_harmonics.get_execution_policy(), 
			chroma.get_execution_policy(),
			impulse_harmonic_frequency.get_execution_policy()
			);

	return synthesize_grains( length, grains_per_second, time_scatter, Function<Second, Audio>( [&]( Second t )
		{ 
		const Audio impulse = synthesize_impulse( impulse_harmonic_frequency( t ), num_harmonics( t ), chroma( t ), sample_rate );
		return impulse.synthesize_grains_repeat(
			trainlet_length( t ),
			freq,
			0,
			trainlet_gain_envelope,
			sample_rate 
			).stereo_spatialize( position( t ) );
		}, ex_pol ) );
	}

Audio Audio::synthesize_granulation( 
	Second length, 
	const Function<Second, float> & grains_per_second, 
	const Function<Second, float> & time_scatter, 
	const Function<Second, Second> & time_selection, 
	const Function<Second, Second> & grain_length,
	const AudioMod & mod
	) const
	{
	if( is_null() ) return Audio::create_null();

	const auto selection_sampled 	= time_selection.sample( 0, time_to_frame( length ), frame_to_time( 1 ) );
	const auto grain_length_sampled = grain_length	.sample( 0, time_to_frame( length ), frame_to_time( 1 ) );

	auto grain_generator = [&]( const Audio & in, float t ) -> Audio
		{
		const Second selection_c = selection_sampled[ time_to_frame( t ) ];
		const Second grain_length_c = grain_length_sampled[ time_to_frame( t ) ];
		Audio grain = in.cut( selection_c, selection_c + grain_length_c );
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

Audio Audio::synthesize_psola(
	Second length, 
	const Function<Second, Second> & time_selection,
	const AudioMod & mod
	) const
	{
	auto freq = get_frequency_envelope();

	//auto grains_per_second = [&]( Second t ){ return 1.0f / freq( t ); };
	auto grain_length = [&]( Second t ){ return 2.0f / freq( t ); };
	AudioMod composition_mod = [&]( Audio & a, Second t )
		{
		if( ! mod.is_null() )
			mod( a, t );

		a.modify_volume_in_place( [&a]( Second t ){ return Windows::hann( t / a.get_length() ); } );
		};

	return synthesize_granulation( 
		length,
		freq,
		0,
		time_selection,
		grain_length,
		composition_mod 
		);
	}

}
