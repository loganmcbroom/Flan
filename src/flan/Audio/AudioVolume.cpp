#include "flan/Audio/Audio.h"

using namespace flan;

Audio Audio::modify_volume( 
	const Function<Second, float> & volume_level 
	) const
	{
	if( is_null() ) return Audio::create_null();
	Audio out = copy();
	out.modify_volume_in_place( volume_level );
	return out;
	}

Audio& Audio::modify_volume_in_place( 
	const Function<Second, Amplitude> & gain 
	)
	{
	auto gain_sampled = sample_function_over_domain( gain );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		flan::for_each_i( get_num_frames(), ExecutionPolicy::Parallel_Unsequenced, [&]( Frame frame )
			{
			get_sample( channel, frame ) *= gain_sampled[frame];
			} );
	return *this;
	}

Audio Audio::set_volume( 
	const Function<Second, Amplitude> & level 
	) const
	{
	if( is_null() ) return Audio::create_null();
	Audio out = copy();
	out.set_volume_in_place( level );
	return out;
	}

Audio& Audio::set_volume_in_place(
	const Function<Second, Amplitude> & level
	)
	{
	if( is_null() ) return *this;

	// Divide by get_max_sample_magnitude to normalize, multiply by level to set
	const Sample max_mag = get_max_sample_magnitude();
	if( max_mag == 0 ) return *this;
	else modify_volume_in_place( [&]( float t ){ return level(t) / max_mag; } );
	}

Audio Audio::fade( 
	Second start, 
	Second end, 
	Interpolator interp 
	) const
	{
	if( is_null() ) return Audio::create_null();
	return fade_frames( time_to_frame( start ), time_to_frame( end ), interp );
	}

Audio& Audio::fade_in_place( 
	Second start, 
	Second end, 
	Interpolator interp 
	)
	{
	if( is_null() ) { std::cout << "Null input" << std::endl; return *this; }
	fade_frames_in_place( time_to_frame( start ), time_to_frame( end ), interp );
	return *this;
	}

Audio Audio::fade_frames( 
	Frame start, 
	Frame end, 
	Interpolator interp 
	) const
	{
	if( is_null() ) return Audio::create_null();
	Audio out = copy();
	out.fade_frames_in_place( start, end, interp );
	return out;
	}

Audio& Audio::fade_frames_in_place( 
	Frame start, 
	Frame end, 
	Interpolator interp 
	)
	{
	if( is_null() ) { std::cout << "Null input" << std::endl; return *this; }

	// Input validation
	start = std::max( 0, start );
	end   = std::max( 0, end   );
	if( start + end > get_num_frames() )
		{
		const float scale = float( get_num_frames() ) / ( start + end );
		start = std::floor( start * scale );
		end	= std::floor( end * scale );
		}
	if( start == 0 && end == 0 ) return *this;

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		for( Frame frame = 0; frame < start; ++frame )
			{
			const float fade_amount = interp( float( frame ) / start );
			get_sample( channel, frame ) *= fade_amount;
			}
		for( Frame frame = 0; frame < end; ++frame )
			{
			const float fade_amount = interp( float( frame ) / end );
			get_sample( channel, get_num_frames() - 1 - frame ) *= fade_amount;
			}
		}
	return *this;
	}

Audio Audio::invert_phase(
	) const
	{
	if( is_null() ) return Audio::create_null();
	Audio out = copy();
	std::for_each( FLAN_PAR_UNSEQ out.get_buffer().begin(), out.get_buffer().end(), []( Sample & s ){ s = -s; } );
	return out;
	}

Audio Audio::waveshape( 
	const Function< std::pair<Second, Sample>, Sample > & shaper,
	uint16_t oversample_factor
	) const
	{
	if( is_null() ) return Audio::create_null();

	Audio oversampled = resample( get_sample_rate() * oversample_factor );
	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		runtime_execution_policy_handler( shaper.get_execution_policy(), [&]( auto policy )
			{
			std::for_each( FLAN_POLICY iota_iter( 0 ), iota_iter( oversampled.get_num_frames() ), [&]( Frame frame )
				{ 
				Sample & s = oversampled.get_sample( channel, frame );
				s = shaper( std::pair( oversampled.frame_to_time( frame ), s ) );
				} );
			} );
		}
	return oversampled.resample( get_sample_rate() );
	}

Audio Audio::add_moisture(
	const Function<Second, Amplitude> & amount,
	const Function<Second, Frequency> & frequency,
	const Function<Second, float> & skew,
	const Function<Second, Amplitude> & waveform
	) const
	{
	auto amount_sampled = sample_function_over_domain( amount );
	auto frequency_sampled = sample_function_over_domain( frequency );
	auto skew_sampled = sample_function_over_domain( skew );

	return waveshape( [&]( std::pair<Second, Sample> ts ) -> Sample
		{ 
		const float amount_c 	= amount_sampled	[time_to_frame(ts.first)];
		const float frequency_c = frequency_sampled	[time_to_frame(ts.first)];
		const float skew_c 		= skew_sampled		[time_to_frame(ts.first)];

		const float power = ts.second >= 0 ? std::pow( ts.second, skew_c ) : -std::pow( -ts.second, skew_c );
		return ts.second + amount_c * ts.second * waveform( pi2 * frequency_c * power ); 
		} );
	}

Audio Audio::compress( 
	const Function<Second, Decibel> & threshold, 
	const Function<Second, float> & ratio, 
	const Function<Second, Second> & attack, 
	const Function<Second, Second> & release, 
	const Function<Second, Decibel> & knee_width, 
	const Audio * sidechain_source 
	) const
	{
	/*
	For details on compression, both in general and as referenced in this function:
	"Digital Dynamic Range Compressor Design — A Tutorial and Analysis"
	https://www.eecs.qmul.ac.uk/~josh/documents/2012/GiannoulisMassbergReiss-dynamicrangecompression-JAES2012.pdf
	*/

	if( is_null() ) return Audio::create_null();

	if( sidechain_source == nullptr )
		sidechain_source = this;

	// All channels need to be compressed equally, so the volume control input signal used is the max signal over channels
	std::vector<Sample> channel_max( sidechain_source->get_num_frames(), 0 );
	for( Channel c = 0; c < sidechain_source->get_num_channels(); ++c )
		for( Frame f = 0; f < sidechain_source->get_num_frames(); ++f )
			if( channel_max[f] < sidechain_source->get_sample( c, f ) )
				channel_max[f] = sidechain_source->get_sample( c, f );

	Audio out = copy();

	// Sample all inputs
	auto threshold_sampled 	= threshold .sample( 0, get_num_frames(), frame_to_time( 1 ) );
	auto ratio_sampled 		= ratio     .sample( 0, get_num_frames(), frame_to_time( 1 ) );
	auto attack_sampled 	= attack    .sample( 0, get_num_frames(), frame_to_time( 1 ) );
	auto release_sampled 	= release   .sample( 0, get_num_frames(), frame_to_time( 1 ) );
	auto knee_width_sampled = knee_width.sample( 0, get_num_frames(), frame_to_time( 1 ) );

	// (4)
	auto gain_computer = [&]( Decibel x_G, Decibel threshold, Decibel knee_width, float ratio ) -> Decibel
		{
		const Decibel overshoot = x_G - threshold;
		if( overshoot <= -knee_width/2.0f ) // Before knee
			return x_G;
		else if( overshoot >= knee_width / 2.0f ) // After knee
			return x_G + overshoot * ( 1 / ratio - 1 );
		else // In the knee
			{
			const Decibel z = overshoot + knee_width/2.0f;
			return x_G + ( 1 / ratio - 1 ) * z*z / ( 2.0f * knee_width );
			}
		};

	// (7)
	auto time_to_alpha = [sr = get_sample_rate()]( Second t ){ return std::exp( -1.0f / ( t * sr ) ); };

	// (17)
	// The term "peak detector" is used in the literature. It is a filter.
	auto smooth_decouple_peak_detector = [&]( float & y_1, float & y_L, Frame f, float x_L )
		{
		const float a_A = time_to_alpha( attack_sampled[ f ] );
		const float a_R = time_to_alpha( release_sampled[ f ] );
		y_1 = std::max( x_L, a_R * y_1 + ( 1.0f - a_R ) * x_L );
		y_L = a_A * y_L + ( 1.0f - a_A ) * y_1;
		return y_L;
		};

	flan::for_each_i( get_num_channels(), ExecutionPolicy::Parallel_Unsequenced, [&]( Channel channel)
		{
		// Peak detector buffers
		float y_1 = 0;
		float y_L = 0;
		for( Frame frame = 0; frame < get_num_frames(); ++frame )
			{
			// (23)
			const Sample x = channel_max[frame];
			const Decibel x_G = 20.0f * std::log10( std::max( std::abs( x ), 1e-6f ) ); // Max is to avoid -inf from log
			const Decibel y_G = gain_computer( x_G, 
				threshold_sampled[ frame ], 
				knee_width_sampled[ frame ], 
				ratio_sampled[ frame ] );
			const Decibel x_L = x_G - y_G;
			const Decibel c_dB = -smooth_decouple_peak_detector( y_1, y_L, frame, x_L );

			const Sample c = std::pow( 10.0f, c_dB / 20.0f );
			out.get_sample( channel, frame ) *= c;
			}
		} );

	return out;
	}
