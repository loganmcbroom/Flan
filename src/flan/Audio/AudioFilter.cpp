#include "flan/Audio/Audio.h"

/*
https://ia601900.us.archive.org/5/items/the-art-of-va-filter-design-rev.-2.1.2/VAFilterDesign_2.1.2.pdf#chapter.10
*/

using namespace flan;

using Pole = std::complex<float>;
using Mix_1pole = std::array<float,2>;
using Mix_2pole = std::array<float,3>;
using Mix_Func_1pole = Function<Second, Mix_1pole>;
using Mix_Func_2pole = Function<Second, Mix_2pole>;

//===============================================================================================================================
// Utility
//===============================================================================================================================

Frequency prewarp( Frequency w, float T_half ) 
	{
	/*
	Filter digitization bends the frequency response near nyquist down slightly.
	Rather than have this error present in full at nyquist and not at all near zero, we can increase the desired frequency
	in such a way that after this "bending" occurs, the cutoff position of the filter is right on the money.
	This does create additional response error around 0, but we can only minimize that error around a single frequency,
	and the cutoff is the most audible.
	*/

	return std::tan( T_half*w ) / T_half;
	}

std::vector<Pole> generate_butterworth_type1_poles( uint16_t N )
	{
	// Get type 1 roots with unit cutoff. 
	// Note, we are only listing the roots above the x-axis. For each non-real root, there is an unlisted complex conjugate root.
	std::vector<Pole> poles;
	for( Index i = 0; i < std::floor( N / 2.0f ); ++i )
		{
		const Radian delta = pi2 / ( N * 2 );
		const Radian theta = delta * i + pi/2.0f + delta/2.0f;
		poles.push_back( std::exp( Pole( 0, theta ) ) );
		};
	return poles;
	}



//===============================================================================================================================
// 1-pole base
//===============================================================================================================================

struct Filter_1Pole {
	Filter_1Pole( FrameRate sr )
		: s( 0 )
		// Note the factor of 2pi. The reference book uses a non-standard Forier transform definition, much to my annoyance. This fixes that.
		// Saving a few cycles, the factor of 1/2 is baked into T, taking the 2pi down to pi.
		, T_half( pi / sr )
		{
		}

	std::array<Sample, 2> process_sample( Sample x, Frequency cutoff_unwarped, bool use_prewarp = true )
		{
		// See section 3.10 for the details.
		
		const Frequency w = use_prewarp? prewarp( cutoff_unwarped, T_half ) : cutoff_unwarped;
		
		const float g = w * T_half; 
		const float G = g / ( 1 + g );
		const float v = G * ( x - s );
		const Sample lp = v + s;
		s = lp + v;

		return { lp, x - lp };
		}

	Sample process_sample_and_mix( Sample x, Frequency cutoff_unwarped, Mix_1pole mix, bool use_prewarp = true )
		{
		const auto filtered = process_sample( x, cutoff_unwarped, use_prewarp );
		return filtered[0]*mix[0] + filtered[1]*mix[1];
		}

	Sample s;
	const float T_half;
};

// std::vector<Audio> base_filter_1pole_multimode( 
// 	const Audio & me,
// 	const Function<Second, Frequency> & cutoff,
// 	const std::vector<std::array<float, 2>> & mixers
// 	)
// 	{
// 	std::vector<Audio> outs;
// 	for( auto & _ : mixers ) outs.emplace_back( me.get_format() );

// 	for( Channel channel = 0; channel < me.get_num_channels(); ++channel )
// 		{
// 		Filter_1Pole filter( me.get_sample_rate() );
// 		for( Frame frame = 0; frame < me.get_num_frames(); ++frame )
// 			{
// 			auto filtered = filter.processSample( me.get_sample( channel, frame ), cutoff( me.frame_to_time( frame ) ) );

// 			for( Index i = 0; i < mixers.size(); ++i )
// 				{
// 				const std::array<float, 2> m = mixers[i];
// 				outs[i].get_sample( channel, frame ) = m[0]*filtered[0] + m[1]*filtered[1];
// 				}
// 			}
// 		}

// 	return outs;
// 	}

Audio base_filter_1pole_selector( 
	const Audio & me,
	const Function<Second, Frequency> & cutoff,
	size_t i
	)
	{
	auto cutoff_sampled = me.sample_function_over_domain( cutoff );
	cutoff_sampled.for_each( [&]( auto & c ){ c = std::clamp( c, 1.0f, me.get_sample_rate()/2.0f ); } );

	Audio out( me.get_format() );
	for( Channel channel = 0; channel < me.get_num_channels(); ++channel )
		{
		Filter_1Pole filter( me.get_sample_rate() );
		for( Frame frame = 0; frame < me.get_num_frames(); ++frame )
			out.get_sample( channel, frame ) = filter.process_sample( me.get_sample( channel, frame ), cutoff_sampled[frame] )[i];
		}
	return out;
	}

Audio base_filter_1pole_lowpass( 
	const Audio & me,
	const Function<Second, Frequency> & cutoff 
	)
	{
	return base_filter_1pole_selector( me, cutoff, 0 );
	}

Audio base_filter_1pole_highpass( 
	const Audio & me,
	const Function<Second, Frequency> & cutoff
	)
	{
	return base_filter_1pole_selector( me, cutoff, 1 );
	}



//===============================================================================================================================
// 2-pole base
//===============================================================================================================================

struct Filter_2Pole {
	Filter_2Pole( FrameRate sr )
		: s1( 0 )
		, s2( 0 )
		// Note the factor of 2pi. The reference book uses a non-standard Forier transform definition, much to my annoyance. This fixes that.
		// Saving a few cycles, the factor of 1/2 is baked into T, taking the 2pi down to pi.
		, T_half( pi / sr )
		{
		}

	std::array<Sample, 3> process_sample( Sample x, Frequency cutoff_unwarped, float R )
		{
		//See section 4.4 for the implementation.

		const Frequency w = prewarp( cutoff_unwarped, T_half );

		const float g = w * T_half;
		const float g1 = 2.0f*R + g;
		const float d = 1.0f / ( 1.0f + 2.0f*R*g + g*g );
		const float hp = ( x - g1*s1 - s2 ) * d;
		const float v1 = g*hp;
		const float bp = v1 + s1;
		s1 = bp + v1;
		const float v2 = g*bp;
		const float lp = v2 + s2;
		s2 = lp + v2;

		return { lp, bp*2*R, hp };
		}

	Sample process_sample_and_mix( Sample x, Frequency cutoff_unwarped, float R, Mix_2pole mix )
		{
		const auto filtered = process_sample( x, cutoff_unwarped, R );
		return filtered[0]*mix[0] + filtered[1]*mix[1] + filtered[2]*mix[2];
		}

	Sample s1, s2;
	const float T_half;
};

// This base function is needed because R and w sometimes use one another
// Having a single function return both would make a bad forward interface but saves several trig calls per frame sometimes
Audio base_filter_2pole_selector_single_func(
	const Audio & me,
	const Function<Second, std::pair<Frequency, float>> & cutoff_and_damping,
	size_t i
	)
	{
	const auto cutoff_and_damping_sampled = me.sample_function_over_domain( cutoff_and_damping );

	Audio out( me.get_format() ); 

	for( Channel channel = 0; channel < me.get_num_channels(); ++channel )
		{
		Filter_2Pole filter( me.get_sample_rate() );
		for( Frame frame = 0; frame < me.get_num_frames(); ++frame )
			{
			const auto w_R_c = cutoff_and_damping_sampled[frame];
			out.get_sample( channel, frame ) = filter.process_sample( me.get_sample( channel, frame ), w_R_c.first, w_R_c.second )[i];
			}
		}

	return out;
	}
	
Audio base_filter_2pole_selector_from_pole(
	const Audio & me,
	const Function<Second, Pole> & pole,
	size_t i
	)
	{
	return base_filter_2pole_selector_single_func( me,
		Function<Second, std::pair<Frequency, float>>( [&]( Second t )
			{
			const Pole pole_c = pole( t );
			const Magnitude pole_abs = std::abs( pole_c );
			return std::make_pair( pole_abs, -pole_c.real() / pole_abs ); 
			}, pole.get_execution_policy() ),
		i );
	}

Audio base_filter_2pole_selector(
	const Audio & me,
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping,
	size_t i
	)
	{
	return base_filter_2pole_selector_single_func( me, Function<Second, std::pair<Frequency, float>>( 
			[&]( Second t ){ return std::make_pair( cutoff( t ), damping( t ) ); }, lowest_execution( cutoff, damping )
		), i );
	}

Audio base_filter_2pole_lowpass(
	const Audio & me,
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping
	)
	{
	return std::move( base_filter_2pole_selector( me, cutoff, damping, 0 ) );
	}

Audio base_filter_2pole_bandpass(
	const Audio & me,
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping
	)
	{
	return std::move( base_filter_2pole_selector( me, cutoff, damping, 1 ) );
	}

Audio base_filter_2pole_highpass(
	const Audio & me,
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping
	)
	{
	return std::move( base_filter_2pole_selector( me, cutoff, damping, 2 ) );
	}



//===============================================================================================================================
// 1-pole butterworth
//===============================================================================================================================

Audio filter_1pole_repeat_base(
	const Audio & me,
	const Function<Second, Frequency> & cutoff,
	const uint16_t repeats,
	int index
	)
	{
	if( me.is_null() ) return Audio::create_null();
	
	Audio out( me.get_format() );

	auto cutoff_sampled = me.sample_function_over_domain( cutoff );
	cutoff_sampled.for_each( [&]( auto & c ){ c = std::clamp( c, 1.0f, me.get_sample_rate()/2.0f ); } );

	for_each_i( me.get_num_channels(), ExecutionPolicy::Parallel_Unsequenced, [&]( Channel channel )
		{
		std::vector<Filter_1Pole> filter_1poles( repeats, me.get_sample_rate() );
		for( Frame frame = 0; frame < me.get_num_frames(); ++frame )
			{
			for( Index i = 0; i < repeats; ++i )
				{
				const Audio & source = i == 0? me : out;
				out.get_sample( channel, frame ) = filter_1poles[i].process_sample( source.get_sample( channel, frame ), cutoff_sampled[frame] )[index];
				}
			}
		} );

	return out;
	}

Audio Audio::filter_1pole_repeat_low(
	const Function<Second, Frequency> & cutoff,
	const uint16_t repeats
	) const
	{
	return filter_1pole_repeat_base( *this, cutoff, repeats, 0 );
	}

Audio Audio::filter_1pole_repeat_high(
	const Function<Second, Frequency> & cutoff,
	const uint16_t repeats
	) const
	{
	return filter_1pole_repeat_base( *this, cutoff, repeats, 1 );
	}
	

Audio base_filter_1pole_butterworth_selector(
	const Audio & me,
	uint16_t order,
	const Function<Second, Frequency> & cutoff,
	bool lowpass
	)
	{
	/*
	See section 8.6 for details.
	*/
	if( order == 0 ) return me.copy();
	const bool even_order = order % 2 == 0;

	const std::vector<Pole> poles = generate_butterworth_type1_poles( order );

	auto cutoff_sampled = me.sample_function_over_domain( cutoff );
	cutoff_sampled.for_each( [&]( auto & c ){ c = std::clamp( c, 1.0f, me.get_sample_rate()/2.0f ); } );

	Audio out( me.get_format() );

	for( Channel channel = 0; channel < me.get_num_channels(); ++channel )
		{
		Filter_1Pole filter_1pole( me.get_sample_rate() );
		std::vector<Filter_2Pole> filter_2poles( poles.size(), me.get_sample_rate() ); 
		for( Frame frame = 0; frame < me.get_num_frames(); ++frame )
			{
			const Frequency w = cutoff_sampled[frame];

			// For odd orders there is a pole at -1
			if( !even_order ) 
				out.get_sample( channel, frame ) = filter_1pole.process_sample( me.get_sample( channel, frame ), w )[ lowpass? 0:1];

			for( Index pole_i = 0; pole_i < poles.size(); ++pole_i )
				{
				const float R = -poles[pole_i].real();
				const Audio & source = even_order && pole_i == 0? me : out;
				out.get_sample( channel, frame ) = filter_2poles[pole_i].process_sample( source.get_sample( channel, frame ), w, R )[ lowpass? 0:2];
				}
			}
		}

	return out;	
	}

Audio Audio::filter_1pole_lowpass(
	const Function<Second, Frequency> & cutoff,
	uint16_t order
	) const
	{
	if( is_null() ) return Audio::create_null();
	return base_filter_1pole_butterworth_selector( *this, order, cutoff, true );
	}

Audio Audio::filter_1pole_highpass(
	const Function<Second, Frequency> & cutoff,
	uint16_t order
	) const
	{
	if( is_null() ) return Audio::create_null();
	return base_filter_1pole_butterworth_selector( *this, order, cutoff, false );
	}

std::vector<Audio> Audio::filter_1pole_split(
	const Function<Second, Frequency> & cutoff,
	uint16_t order
	) const
	{
	/*
	This may have you saying "that doesn't seem quite right".
	I am tired of making filters and I'm choosing to not write a perfect L-R crossover filter.
	Lets just agree this works pretty well and move on.
	*/

	// Classic janky forced single call for functional input
	auto cutoff_sampled = sample_function_over_domain( cutoff );
	cutoff_sampled.for_each( [&]( auto & c ){ c = std::clamp( c, 1.0f, get_sample_rate()/2.0f ); } );
	
	auto cutoff_no_resample_func = [&]( Second t ){ return cutoff_sampled[time_to_frame( t )]; };

	if( order <= 1 )
		{
		std::vector<Audio> outs; 
		outs.push_back( filter_1pole_lowpass ( cutoff_no_resample_func, 1 ) ); 
		outs.push_back( filter_1pole_highpass( cutoff_no_resample_func, 1 ) ); 
		return outs;
		}

	std::vector<Audio> outs;
	outs.push_back( 
		 filter_1pole_lowpass( cutoff_no_resample_func, order ) 
		.filter_1pole_lowpass( cutoff_no_resample_func, order ) 
		);
	outs.push_back( 
		 filter_1pole_highpass( cutoff_no_resample_func, order ) 
		.filter_1pole_highpass( cutoff_no_resample_func, order ) 
		);

	return outs;	
	}

//===============================================================================================================================
// 1-pole butterworth shelving
//===============================================================================================================================

Audio base_filter_1pole_butterworth_tilt(
	const Audio & me,
	uint16_t order,
	const Function<Second, Frequency> & cutoff,
	const Function<Second, Decibel> & gain
	)
	{
	/*
	See section 10.2 for the background for the shelving filters.
	For a 1-pole transfer function G_w(s)=w/(s+w), simplifying G(s/M)/G(MS) gives G_Mw(s) + M^2 * ( s / (s + w) ).
	That is, a cutoff scaled copy of a lowpass and an amplitude scaled copy of the corresponding highpass.
	This produces a high shelving filter which can then be scaled to obtain a low shelf.

	Section 10.4 covers 2-pole shelving. With both those in hand, arbitrary shelving filters can be constructed as cascades.
	*/

	if( order == 0 ) return me.copy();
	const bool even_order = order % 2 == 0;

	const std::vector<Pole> poles = generate_butterworth_type1_poles( order );

	auto cutoff_sampled = me.sample_function_over_domain( cutoff );
	cutoff_sampled.for_each( [&]( auto & c ){ c = std::clamp( c, 1.0f, me.get_sample_rate()/2.0f ); } );
	const auto gain_sampled = me.sample_function_over_domain( gain );

	Audio out( me.get_format() );

	for( Channel channel = 0; channel < me.get_num_channels(); ++channel )
		{
		Filter_1Pole filter_1pole( me.get_sample_rate() );
		std::vector<Filter_2Pole> filter_2poles( poles.size(), me.get_sample_rate() ); 
		for( Frame frame = 0; frame < me.get_num_frames(); ++frame )
			{
			//const Second t = me.frame_to_time( frame );
			// Division by 2*order allows specifying the filter in terms of total shelf gain rather than gain per octave.
			const float M = decibel_to_amplitude( gain_sampled[frame] / ( 2 * order ) );
			const float M2 = M * M;
			const Frequency w = M * cutoff_sampled[frame];

			// For odd orders there is a pole at -1
			if( !even_order ) 
				{
				const auto filtered = filter_1pole.process_sample( me.get_sample( channel, frame ), w );
				out.get_sample( channel, frame ) = filtered[0] * M + filtered[1] / M;
				}

			for( Index pole_i = 0; pole_i < poles.size(); ++pole_i )
				{
				const float R = poles[pole_i].real() / w;
				const Audio & source = even_order && pole_i == 0? me : out;
				const auto filtered = filter_2poles[pole_i].process_sample( source.get_sample( channel, frame ), w, R );
				out.get_sample( channel, frame ) = filtered[0] / M2 + filtered[1] + filtered[2] * M2;
				}
			}
		}

	return out;
	}

Audio Audio::filter_1pole_lowshelf(
	const Function<Second, Frequency> & cutoff,
	const Function<Second, Decibel> & gain,
	uint16_t order
	) const
	{
	if( is_null() ) return Audio::create_null();
	Audio tilt = base_filter_1pole_butterworth_tilt( *this, order, cutoff, gain );
	tilt.modify_volume_in_place( Function<float, float>( [&]( Second t ){ return decibel_to_amplitude( gain( t ) / 2 ); }, gain.get_execution_policy() ) );
	return tilt;
	}

Audio Audio::filter_1pole_highshelf(
	const Function<Second, Frequency> & cutoff,
	const Function<Second, Decibel> & gain,
	uint16_t order
	) const
	{
	if( is_null() ) return Audio::create_null();
	Audio tilt = base_filter_1pole_butterworth_tilt( *this, order, cutoff, Function<float, float>( [&]( Second t ){ return -gain( t ); }, gain.get_execution_policy() ) );
	tilt.modify_volume_in_place( Function<float, float>( [&]( Second t ){ return decibel_to_amplitude( gain( t ) / 2 ); }, gain.get_execution_policy() ) );
	return tilt;
	}



//===============================================================================================================================
// 2-pole butterworth
//===============================================================================================================================

Audio base_filter_2pole_butterworth_selector(
	const Audio & me,
	uint16_t order,
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping, 
	size_t i
	)
	{
	if( order == 0 ) return me.copy();
	const bool even_order = order % 2 == 0;

	const std::vector<Pole> poles = generate_butterworth_type1_poles( order );

	auto cutoff_sampled = me.sample_function_over_domain( cutoff );
	cutoff_sampled.for_each( [&]( auto & c ){ c = std::clamp( c, 1.0f, me.get_sample_rate()/2.0f ); } );
	const auto damping_sampled = me.sample_function_over_domain( damping );

	Audio out( me.get_format() );

	for( Channel channel = 0; channel < me.get_num_channels(); ++channel )
		{
		Filter_2Pole filter_1pole( me.get_sample_rate() );
		std::vector<std::array<Filter_2Pole, 2>> filter_2poles( poles.size(), { me.get_sample_rate(), me.get_sample_rate() } ); 
		for( Frame frame = 0; frame < me.get_num_frames(); ++frame )
			{
			const Frequency w = cutoff_sampled[frame];
			const float R = damping_sampled[frame];
			const Radian alpha = std::acos( R ) / order;
			
			// For odd orders there is a pole at -1 that will be split into two reciprocal poles around 0.
			// These poles can be handled with a single svf.
			if( !even_order ) 
				{
				const float real_pole_R = std::cos( alpha );
				out.get_sample( channel, frame ) = filter_1pole.process_sample( me.get_sample( channel, frame ), w, real_pole_R )[i];
				}

			// Each complex pole generated thus far will be split into two poles, each of which is then representing also its conjugate.
			// In total, each pole in runs two 2-pole applications
			for( Index pole_i = 0; pole_i < poles.size(); ++pole_i )
				{
				const std::complex<float> pole_scaler = R > 1 ? 
					std::pow( R + std::sqrt( R*R - 1.0f ), 1.0f / order ) : 
					std::exp( std::complex<float>( 0, -alpha ) );

				const Pole p_w = poles[pole_i] * w;
				const Audio & source = even_order && pole_i == 0? me : out;

				const Pole p1 = p_w * pole_scaler;
				const Frequency p1_w = std::abs( p1 );
				const float p1_R = -p1.real() / p1_w;
				out.get_sample( channel, frame ) = filter_2poles[pole_i][0].process_sample( source.get_sample( channel, frame ), p1_w, p1_R )[i];

				const Pole p2 = p_w / pole_scaler;
				const Frequency p2_w = std::abs( p2 );
				const float p2_R = -p2.real() / p2_w;
				out.get_sample( channel, frame ) = filter_2poles[pole_i][1].process_sample( out   .get_sample( channel, frame ), p2_w, p2_R )[i];
				}
			}
		}

	return out;	
	}

Audio Audio::filter_2pole_lowpass(
	const Function<Second, Frequency> & w,
	const Function<Second, float> & R,
	uint16_t N
	) const
	{
	if( is_null() ) return Audio::create_null();
	return base_filter_2pole_butterworth_selector( *this, N, w, R, 0 );
	}

Audio Audio::filter_2pole_bandpass(
	const Function<Second, Frequency> & w,
	const Function<Second, float> & R,
	uint16_t N
	) const
	{
	if( is_null() ) return Audio::create_null();
	return base_filter_2pole_butterworth_selector( *this, N, w, R, 1 );
	}

Audio Audio::filter_2pole_highpass(
	const Function<Second, Frequency> & w,
	const Function<Second, float> & R,
	uint16_t N
	) const
	{
	if( is_null() ) return Audio::create_null();
	return base_filter_2pole_butterworth_selector( *this, N, w, R, 2 );
	}

Audio Audio::filter_2pole_notch(
	const Function<Second, Frequency> & w,
	const Function<Second, float> & R,
	uint16_t N
	) const
	{
	if( is_null() ) return Audio::create_null();
	Audio bp = filter_2pole_bandpass( w, R, N );
	bp.modify_volume_in_place( -1 );
	return Audio::mix( bp, *this );
	}


//===============================================================================================================================
// 2-pole butterworth shelving
//===============================================================================================================================

Audio base_filter_2pole_butterworth_tilt(
	const Audio & me,
	uint16_t order,
	const std::function<Frequency (Second, float)> & cutoff,
	const std::function<float (Second, float)> & damping,
	const Function<Second, Decibel> & gain,
	std::function<Sample ( Mix_2pole, float )> mix_filtered_sample
	)
	{
	if( order == 0 ) return me.copy();
	const bool even_order = order % 2 == 0;

	const std::vector<Pole> poles = generate_butterworth_type1_poles( order );

	// Sampling inputs, some requiring manual sampling due to the unusual function signature
	const FunctionSample<Decibel> gain_sampled = me.sample_function_over_domain( gain );
	const FunctionSample<Amplitude> Ms = gain_sampled.transform( [&]( Decibel dB ){ return decibel_to_amplitude( dB / ( 2 * order ) ); } );
	std::vector<Frequency> cutoff_sampled( me.get_num_frames() );
	std::vector<float> damping_sampled( me.get_num_frames() );
	for( Frame frame = 0; frame < me.get_num_frames(); ++frame )
		{
		cutoff_sampled[frame] = cutoff( me.frame_to_time( frame ), Ms[frame] );
		damping_sampled[frame] = damping( me.frame_to_time( frame ), Ms[frame] );
		}
	// Note, mix_filtered_sample doesn't need to be sampled as it's not a part of the public filter interface.

	Audio out( me.get_format() );

	for( Channel channel = 0; channel < me.get_num_channels(); ++channel )
		{
		Filter_2Pole filter_1pole( me.get_sample_rate() );
		std::vector<std::array<Filter_2Pole, 2>> filter_2poles( poles.size(), { me.get_sample_rate(), me.get_sample_rate() } ); 
		for( Frame frame = 0; frame < me.get_num_frames(); ++frame )
			{
			const float M = Ms[frame];
			const float M2 = M * M;
			const Frequency w = cutoff_sampled[frame];
			const float R = damping_sampled[frame];
			const Radian alpha = std::acos( R ) / order;
			
			// For odd orders there is a pole at -1 that will be split into two reciprocal poles around 0.
			// These poles can be handled with a single svf.
			if( !even_order ) 
				{
				const float real_pole_R = std::cos( alpha );
				const auto filtered = filter_1pole.process_sample( me.get_sample( channel, frame ), w, real_pole_R );
				out.get_sample( channel, frame ) = mix_filtered_sample( filtered, M2 );
				}

			// Each complex pole generated thus far will be split into two poles, each of which is then representing also its conjugate.
			// In total, each pole in runs two 2-pole applications
			for( Index pole_i = 0; pole_i < poles.size(); ++pole_i )
				{
				const std::complex<float> pole_scaler = R > 1 ? 
					std::pow( R + std::sqrt( R*R - 1.0f ), 1.0f / order ) : 
					std::exp( std::complex<float>( 0, -alpha ) );

				const Pole p_w = poles[pole_i] * w;
				const Audio & source = even_order && pole_i == 0? me : out;

				const Pole p1 = p_w * pole_scaler;
				const Frequency p1_w = std::abs( p1 );
				const float p1_R = -p1.real() / p1_w;
				const auto filtered_1 = filter_2poles[pole_i][0].process_sample( source.get_sample( channel, frame ), p1_w, p1_R );
				out.get_sample( channel, frame ) = mix_filtered_sample( filtered_1, M2 );

				const Pole p2 = p_w / pole_scaler;
				const Frequency p2_w = std::abs( p2 );
				const float p2_R = -p2.real() / p2_w;
				const auto filtered_2 = filter_2poles[pole_i][1].process_sample( out   .get_sample( channel, frame ), p2_w, p2_R );
				out.get_sample( channel, frame ) = mix_filtered_sample( filtered_2, M2 );
				}
			}
		}

	return out;	
	}

Audio Audio::filter_2pole_lowshelf(
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping,
	const Function<Second, Decibel> & gain,
	uint16_t order
	) const
	{
	if( is_null() ) return Audio::create_null();
	Audio tilt = base_filter_2pole_butterworth_tilt( *this, order, 
		[&]( Second t, float M ){ return cutoff( t ) * M; }, 
		[&]( Second t, float ){ return damping( t ); }, 
		Function<float, float>( [&]( Second t ){ return gain( t ) / 2.0f; }, gain.get_execution_policy() ), 
		[&]( Mix_2pole f, float M2 ){ return f[0] / (M2*M2) + f[1] / M2 + f[2]; } );
	//tilt.modify_volume_in_place( [&]( Second t ){ return decibel_to_amplitude( -half_gain( t ) ); } );
	return tilt;
	}

Audio Audio::filter_2pole_bandshelf(
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping,
	const Function<Second, Decibel> & gain,
	uint16_t order
	) const
	{
	if( is_null() ) return Audio::create_null();
	Audio tilt = base_filter_2pole_butterworth_tilt( *this, order, 
		[&]( Second t, float ){ return cutoff( t ); }, 
		[&]( Second t, float M ){ return damping( t ) * M; }, 
		Function<float, float>( [&]( Second t ){ return -gain( t ); }, gain.get_execution_policy() ), 
		[&]( Mix_2pole f, float M2 ){ return f[0] + f[1] / M2 + f[2]; } );
	return tilt;
	}

Audio Audio::filter_2pole_highshelf(
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping,
	const Function<Second, Decibel> & gain,
	uint16_t order
	) const
	{
	if( is_null() ) return Audio::create_null();
	auto half_gain = [&]( Second t ){ return gain( t ) / 2.0f; };
	Audio tilt = base_filter_2pole_butterworth_tilt( *this, order, 
		[&]( Second t, float M ){ return cutoff( t ) * M; }, 
		[&]( Second t, float ){ return damping( t ); }, 
		Function<float, float>( [&]( Second t ){ return gain( t ) / 2.0f; }, gain.get_execution_policy() ), 
		[&]( Mix_2pole f, float M2 ){ return f[0] + f[1] * M2 + f[2] * M2*M2; } );
	//tilt.modify_volume_in_place( [&]( Second t ){ return decibel_to_amplitude( half_gain( t ) ); } );
	return tilt;
	}

//===============================================================================================================================
// Crossover
//===============================================================================================================================

// std::array<Audio, 2> Audio::filter_crossover(
// 	const Function<Second, Frequency> & cutoff,
// 	uint16_t order
// 	) const
// 	{
// 	// See section 10.8

// 	std::array<Audio, 2> out = { get_format(), get_format() };

// 	return out;
// 	}


//===============================================================================================================================
// Multinotch 
//===============================================================================================================================

// Audio filter_allpass(
// 	const Audio & me,
// 	std::function< Sample ( Sample, Second ) > & allpass,
// 	const std::functon< float ( Second ) > & feedback,
// 	bool invert
// 	)
// 	{
// 	Audio allpass = me.get_format();

// 	for( Channel channel = 0; channel < me.get_num_channels(); ++channel )
// 		{
// 		for( Frame frame = 0; frame < me.get_num_frames(); ++frame )
// 			{
// 			const Sample x = me.get_sample( channel, frame );
// 			const Sample G = allpass( x, me.frame_to_time( frame ) );	
// 			const Sample H =
// 			allpass.get_sample( channel, frame ) = allpassed;
// 			}
// 		}
// 	}

Audio Audio::filter_1pole_allpass_n(
	const Function<Second, Frequency> & cutoff,
	uint16_t order
	) const
	{
	auto cutoff_sampled = sample_function_over_domain( cutoff );
	cutoff_sampled.for_each( [&]( auto & c ){ c = std::clamp( c, 1.0f, get_sample_rate()/2.0f ); } );

	Audio allpass = Audio::create_from_format( get_format() );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		std::vector<Filter_1Pole> filters( order, get_sample_rate() );
		for( Frame frame = 0; frame < get_num_frames(); ++frame )
			{
			Sample allpassed;
			const Frequency w = cutoff_sampled[frame];
			for( Index filter_i = 0; filter_i < filters.size(); ++filter_i )
				{
				const Sample input_sample = filter_i == 0? get_sample( channel, frame ) : allpassed;
				allpassed = filters[filter_i].process_sample_and_mix( input_sample, w, { 1.0f, -1.0f } );
				}
			allpass.get_sample( channel, frame ) = allpassed;
			}
		}

	return allpass;
	}

Audio Audio::filter_1pole_multinotch(
	uint16_t order,
	const Function<Second, Frequency> & cutoff,
	bool invert,
	const Function<Second, float> & wet_dry
	) const
	{
	if( is_null() ) return Audio::create_null();
	
	// See section 11.1
	Audio allpass = filter_1pole_allpass_n( cutoff, order );
	if( invert ) allpass.modify_volume_in_place( -1 );
	Function<Second, Amplitude> wet_dry_inv( [&]( Second t ){ return 1.0f - wet_dry( t ); }, wet_dry.get_execution_policy() );
	return Audio::mix( { this, &allpass }, {}, { &wet_dry, &wet_dry_inv } );
	}

Audio filter_2pole_allpass_n(
	const Audio & me,
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping,
	uint16_t order
	)
	{
	auto cutoff_sampled = me.sample_function_over_domain( cutoff );
	cutoff_sampled.for_each( [&]( auto & c ){ c = std::clamp( c, 1.0f, me.get_sample_rate()/2.0f ); } );
	const auto damping_sampled = me.sample_function_over_domain( damping );

	Audio allpass = Audio::create_from_format( me.get_format() );

	for( Channel channel = 0; channel < me.get_num_channels(); ++channel )
		{
		std::vector<Filter_2Pole> filters( order, me.get_sample_rate() );
		for( Frame frame = 0; frame < me.get_num_frames(); ++frame )
			{
			Sample allpassed;
			const Frequency w = cutoff_sampled[frame];
			const float R = damping_sampled[frame];
			for( Index filter_i = 0; filter_i < filters.size(); ++filter_i )
				{
				const Sample input_sample = filter_i == 0? me.get_sample( channel, frame ) : allpassed;
				allpassed = filters[filter_i].process_sample_and_mix( input_sample, w, R, { 1.0f, -1.0f, 1.0f } );
				}
			allpass.get_sample( channel, frame ) = allpassed;
			}
		}

	return allpass;
	}

Audio Audio::filter_2pole_multinotch(
	uint16_t order,
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping,
	bool invert,
	const Function<Second, float> & wet_dry
	) const
	{
	if( is_null() ) return Audio::create_null();

	// See section 11.2
	Audio allpass = filter_2pole_allpass_n( *this, cutoff, damping, order );
	if( invert ) allpass.modify_volume_in_place( -1 );
	Function<Second, Amplitude> wet_dry_inv( [&]( Second t ){ return 1.0f - wet_dry( t ); }, wet_dry.get_execution_policy() );
	return Audio::mix( { this, &allpass }, {}, { &wet_dry, &wet_dry_inv } );
	}

Audio Audio::filter_comb(
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & feedback,
	const Function<Second, float> & wet_dry,
	bool invert
	) const
	{
	if( is_null() ) return Audio::create_null();
	
	// See sections 11.5 and 11.6
	// Let k be the feedback amount. Let f by plus or minus one, depending on the comb filter being inverted.
	// Applying 11.6 to a comb filter and letting u_n be the signal entering the comb filter delay element which has a delay of t, we have:
	// (1) u_n = x_n + kfu_(n-t)
	// (2) y_n = 1/2[ u_n + fu_(n-t) ]
	// Using (1) we can compute u_t into an extra buffer
	// Then using (2) we add two copies of the secondary buffer and scale to obtain y_n
	// Note, the code additionally accounts for non-50/50 mixing and inversion, and so is slightly off from the equations.

	const int f = invert? -1 : 1;

	auto cutoff_sampled = sample_function_over_domain( cutoff );
	cutoff_sampled.for_each( [&]( auto & c ){ c = std::clamp( c, 1.0f, get_sample_rate()/2.0f ); } );
	const auto wet_dry_sampled = sample_function_over_domain( wet_dry );
	const auto feedback_sampled = sample_function_over_domain( feedback );

	Audio out = Audio::create_from_format( get_format() );

	std::vector<Sample> u( out.get_num_frames() );
	auto safe_u_access = [&]( Frame frame ) -> Sample
		{ 
		if( frame < 0 || out.get_num_frames() <= frame ) return 0.0f; 
		return u[frame];
		};

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		std::fill( u.begin(), u.end(), 0.0f );
		for( Frame frame = 0; frame < get_num_frames(); ++frame )
			{
			const Frequency w = cutoff_sampled[frame];
			const float k = feedback_sampled[frame];
			const float a = wet_dry_sampled[frame];
			
			const Second delay = 1.0f / (2.0f * w );
			const Sample u_nmt = safe_u_access( frame - time_to_frame( delay ) );

			// (1) u_n = x_n + ku_(n-t)
			u[frame] = get_sample( channel, frame ) + k*f*u_nmt;

			// (2) y_n = 1/2[ u_n + u_(n-t) ]
			out.set_sample( channel, frame, a*u[frame] + (1.0f - a)*f*u_nmt );
			}
		}

	return out;
	}

Audio Audio::filter_1pole_multi_allpass(
	const std::vector<Function<Second, Frequency>> & cutoffs
	) const
	{
	if( cutoffs.empty() ) return Audio::create_null();

	std::vector<FunctionSample<Frequency>> cutoffs_sampled;
	for( int i = 0; i < cutoffs.size(); ++i )
		cutoffs_sampled.push_back( sample_function_over_domain( cutoffs[i] ) );

	Audio allpass = Audio::create_from_format( get_format() );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		std::vector<Filter_1Pole> filters( cutoffs.size(), get_sample_rate() );
		for( Frame frame = 0; frame < get_num_frames(); ++frame )
			{
			Sample allpassed;
			for( Index filter_i = 0; filter_i < filters.size(); ++filter_i )
				{
				const Sample input_sample = filter_i == 0? get_sample( channel, frame ) : allpassed;
				allpassed = filters[filter_i].process_sample_and_mix( input_sample, cutoffs_sampled[filter_i][frame], { 1.0f, -1.0f }, false );
				}
			allpass.get_sample( channel, frame ) = allpassed;
			}
		}

	return allpass;
	}

Audio filter_1pole_multi_allpass_constant(
	const Audio & me,
	const std::vector<Frequency> & cutoffs
	)
	{
	std::vector<Function<Second, Frequency>> cutoff_funcs;
	for( int i = 0; i < cutoffs.size(); ++i )
		cutoff_funcs.push_back( cutoffs[i] );
	return me.filter_1pole_multi_allpass( cutoff_funcs );
	}

// #include <fftw3.h>
// void graph_complex_signal_spectrum( const Audio & real, const Audio & complex )
// 	{
// 	constexpr int N = 262144*2;
// 
//     fftwf_complex * in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
//     fftwf_complex * out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
// 	fftwf_plan p = fftwf_plan_dft_1d( N, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
// 
// 	for( Frame frame = 0; frame < N; ++frame )
// 		{
// 		in[frame][0] = 0;
// 		in[frame][1] = 0;
// 		}
// 	for( Frame frame = 0; frame < real.get_num_frames(); ++frame )
// 		{ 
// 		in[frame][0] = real.get_sample( 0, frame );
// 		in[frame][1] = complex.get_sample( 0, frame );
// 		}
// 
// 	fftwf_execute(p);
// 
// 	Audio spectrum = Audio::create_empty_with_frames( N );
// 	for( Frame frame = 0; frame < spectrum.get_num_frames(); ++frame )
// 		spectrum.get_sample( 0, frame ) = std::log( std::sqrt( out[frame][0]*out[frame][0] + out[frame][1]*out[frame][1] ) + 1 );
// 
// 	auto b = spectrum.set_volume(1).convert_to_graph();
// 	b.save_image( "Temp.bmp" );
// 	system( (std::string("start ") + "Temp.bmp").c_str() );
// 
//     fftwf_destroy_plan(p);
//     fftwf_free(in); 
// 	fftwf_free(out);
// 	}

static std::pair<std::vector<float>, std::vector<float>> phase_diff_network_pole_design( int num_poles, Frequency lower, Frequency upper )
	{
	// See "THE DESIGN OF WIDEBAND ANALOG 90° PHASE DIFFERENCING NETWORKS WITHOUT A LARGE SPREAD OF CAPACITOR VALUES"
	// http://electronotes.netfirms.com/EN168-90degreePDN.PDF

	const Frequency f_l = lower;
	const Frequency f_u = upper;
	const double B = f_u/f_l;
	const double k = std::sqrt(1.0-1.0/(B*B));
	const double L = 0.5*(1.0-std::sqrt(k))/(1.0+std::sqrt(k));
	const double A_p = L + 2.0*std::pow(L,5.0) + 15.0*std::pow(L,9.0);
	const double A = std::exp(std::_Pi*std::_Pi/std::log(A_p));
	//const double E = 1;
	const double n = num_poles;//std::ceil( std::log(E*std::_Pi/720.0)/std::log(A) );

	std::vector<double> phi;
	for( int r = 1; r <= n; ++r )
		phi.push_back( std::_Pi/4.0/n*(2*r-1) );

	std::vector<double> phi_p;
	for( int r = 0; r < phi.size(); ++r)
		{
		const double numerator = (A*A-A*A*A*A*A*A)*std::sin( 4.0*phi[r] );
		const double denominator = 1.0 + (A*A + A*A*A*A*A*A) * std::cos( 4.0*phi[r] );
		phi_p.push_back( std::atan(numerator/denominator) );
		}

	std::vector<double> p;
	for( int r = 0; r < phi.size(); ++r)
		p.push_back( std::sqrt(B) * std::tan(phi[r] - phi_p[r]) * 2.0 * std::_Pi * f_l );

	std::vector<float> p_a;
	std::vector<float> p_b;

	for( int r = 0; r < p.size(); ++r )	
		{
		if( r % 2 == 0 ) p_a.push_back( p[r] );
		else p_b.push_back( p[r] );
		}

	return std::make_pair( p_b, p_a ); // They are backwards, I don't know why but it doesn't matter and it erases the positive spectrum the other way.
	}

std::pair<Audio, Audio> hilbert_phase_diff_network( const Audio & me )
	{
	// Hilbert transform approximation via phase difference network
	// See https://nathan.ho.name/posts/frequency-shifter/
	const auto poles = phase_diff_network_pole_design( 20, 5, 22000 );
	return std::make_pair( 
		filter_1pole_multi_allpass_constant( me, poles.first ),
	 	filter_1pole_multi_allpass_constant( me, poles.second ) );
	}

Audio Audio::halfband_modulate(
	const Function<Second, std::complex<float>> & modulator
	) const
	{
	auto modulator_sampled = sample_function_over_domain( modulator );

	const auto hilberted = hilbert_phase_diff_network( *this );

	Audio out( get_format() );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		flan::for_each_i( hilberted.first.get_num_frames(), ExecutionPolicy::Parallel_Unsequenced, [&]( Frame frame )
			{
			const std::complex<float> mod_c = modulator_sampled[frame];
			const float a = hilberted.first.get_sample( channel, frame ) * mod_c.real();
			const float b = hilberted.second.get_sample( channel, frame ) * mod_c.imag();
			out.get_sample( channel, frame ) = a - b;
			} );
		}

	return out;
	}

Audio Audio::shift_frequency(
	const Function<Second, Frequency> & shift,
	const Frequency low_cutoff
	) const
	{
	const Frequency high_cutoff = get_sample_rate()/2 - 100; // Using exactly nyquist causes feedback in high order filter.

	auto shift_sampled = sample_function_over_domain( shift );

	// Frequencies with 20-20k may be moved supernyquist or sub-dc, which we would not like. This pre-filters those signals away.
	const Audio antialiased = 
		filter_1pole_lowpass( [&]( Second t )
			{ 
			if( shift_sampled[time_to_frame(t)] > 0 )
				return high_cutoff - shift_sampled[time_to_frame(t)];
			else
				return high_cutoff;
			}, 8 )
		.filter_1pole_highpass( [&]( Second t )
			{ 
			if( shift_sampled[time_to_frame(t)] < 0 )
				return low_cutoff-shift_sampled[time_to_frame(t)];
			else
				return low_cutoff;
			}, 8 );

	return antialiased.halfband_modulate( [&]( Second t ){ return std::exp( std::complex<float>( 0.0f, pi2 * shift(t) * t ) ); } );
	}

Audio Audio::halfband_multiply(
	const Audio & modulator
	) const
	{
	auto bandpass_antialias = []( const Audio & a )
		{
		const Frequency low_cutoff = 30;
		const Frequency high_cutoff = a.get_sample_rate()/2 - 2000;
		return a.filter_1pole_lowpass( high_cutoff, 8 ).filter_1pole_highpass( low_cutoff, 8 );
		};

	const auto hilbert1 = hilbert_phase_diff_network( bandpass_antialias( *this ) );
	const auto hilbert2 = hilbert_phase_diff_network( bandpass_antialias( modulator ) );

	Format format = get_format();
	format.num_channels = std::min( get_num_channels(), modulator.get_num_channels() );
	format.num_frames = std::min( get_num_frames(), modulator.get_num_frames() );
	Audio out( format );

	for( Channel channel = 0; channel < out.get_num_channels(); ++channel )
		{
		flan::for_each_i( out.get_num_frames(), ExecutionPolicy::Parallel_Unsequenced, [&]( Frame frame )
			{
			const float a_x = hilbert1.first.get_sample( channel, frame );
			const float a_y = hilbert1.second.get_sample( channel, frame );
			const float b_x = hilbert2.first.get_sample( channel, frame );
			const float b_y = hilbert2.second.get_sample( channel, frame );
			out.get_sample( channel, frame ) = a_x * b_x - a_y * b_y;
			} );
		}
	
	return out;
	}

Audio Audio::phaser_step(
	const Function<Second, float> & lfo_central_cutoff,
	const Function<Second, float> & lfo_amplitude,
	const Function<Second, float> & lfo_frequency,
	float resonance,
	int filter_order,
	Radian initial_phase,
	bool pitch_domain
	) const
	{	
	filter_order = std::clamp( filter_order, 1, 128 );

	auto lfo_central_cutoff_sampled = sample_function_over_domain( lfo_central_cutoff );
	auto lfo_amplitude_sampled = sample_function_over_domain( lfo_amplitude );
	auto lfo_frequency_sampled = sample_function_over_domain( lfo_frequency ); // In cycles per second
	lfo_frequency_sampled.for_each( [&]( float & x ){ x = x / get_sample_rate() * pi2; } ); // Convert to radians per frame
	const std::vector<Radian> phase = lfo_frequency_sampled.exclusive_scan( initial_phase, []( float a, float b ){ return a + b; } );

	return filter_2pole_notch(
		[&]( Second t )
			{ 
			const Frame frame = time_to_frame( t );
			const float cutoff = lfo_central_cutoff_sampled[frame] + lfo_amplitude_sampled[frame] * std::sin( phase[frame] );
			if( pitch_domain ) return std::pow( 2.0f, cutoff ); 
			else return cutoff;
			},
		resonance,
		filter_order );
	}

Audio Audio::easy_phaser(
	const Function<Second, float> & speed,
	const Function<Second, float> & wow_amount,
	int how_many_notches,
	float resonance,
	int filter_order
	) const
	{
	std::random_device rd;
	std::mt19937 rng( rd() );
	std::uniform_real_distribution<> unit_dist( 0.0f, 1.0f );

	filter_order = std::clamp( filter_order, 1, 10 );

	auto speed_sampled = sample_function_over_domain( speed );
	speed_sampled.for_each( [&]( float & x ){ x = x / get_sample_rate() * pi2;} ); // Convert cycles/second to radian/frame
	const std::vector<Radian> phase = speed_sampled.exclusive_scan( 0.0f, [&]( float a, float b ){ return a + b; } ); // Integrate speed to phase
	auto wow_amount_sampled = sample_function_over_domain( wow_amount );

	Audio out = copy();	
	for( int notch = 0; notch < how_many_notches; ++notch )
		{ 
		const float notch_initial_pitch = 9 + notch * 5.0f / how_many_notches + .3 * unit_dist( rng );
		const Radian notch_initial_phase = unit_dist( rng ) * pi2;
		
		const float speed_mod = std::pow( 2.0f, unit_dist( rng ) * 2.0f - 1.0f );

		out = out.filter_2pole_notch(
			[&]( Second t )
				{ 
				const Frame frame = time_to_frame( t );
				const float current_pitch = notch_initial_pitch + wow_amount_sampled[frame] * std::sin( notch_initial_phase + phase[frame]*speed_mod );
				return std::pow( 2.0f, current_pitch ); 
				},
			resonance,
			filter_order );
		}

	return out;
	}

Audio Audio::phaser_texture(
	const Function<Second, float> & effects_per_second, 
	const Function<Second, float> & time_scatter,
	const Function<Second, Second> & effects_length,

	const Function<Second, float> & pitch_start, 
	const Function<Second, float> & pitch_movement,
	const Function<Second, float> & resonance,
	const Function<Second, int> order,

	const Function<Second, float> & dry_wet,
	Second fade_time,
	const Interpolator & interp
	) const
	{
	const ExecutionPolicy lowest_exec = lowest_execution( { 
		pitch_start.get_execution_policy(), 
		pitch_movement.get_execution_policy(),
		resonance.get_execution_policy(),
		order.get_execution_policy() 
		} );

	const Audio wet = texture_effect( effects_per_second, time_scatter, effects_length, 
		AudioMod( [&]( Audio & a, Second t )
			{
			const float start_pitch_c = pitch_start(t);
			const float pitch_movement_c = pitch_movement(t);
			const float R = resonance(t);
			const int order_c = std::clamp( order(t), 1, 128 );
			a = a.filter_2pole_notch( [&]( Second t )
				{ 
				return std::pow( 2.0f, start_pitch_c + pitch_movement_c * interp( t / a.get_length() ) ); 
				}, R, order_c );
			}, lowest_exec ), 
		fade_time );

	const Function<Second, Amplitude> wet_dry( [&]( Second t ){ return 1.0f - dry_wet(t); }, dry_wet.get_execution_policy() );
	const std::vector< const Function<Second, Amplitude> *> gains = { &wet_dry, &dry_wet };
	return Audio::mix( { this, &wet }, {}, gains );
	}