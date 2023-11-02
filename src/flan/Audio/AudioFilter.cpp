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

	std::array<Sample, 2> process_sample( Sample x, Frequency cutoff_unwarped )
		{
		// See section 3.10 for the details.
		
		const Frequency w = prewarp( cutoff_unwarped, T_half );
		
		const float g = w * T_half; 
		const float G = g / ( 1 + g );
		const float v = G * ( x - s );
		const Sample lp = v + s;
		s = lp + v;

		return { lp, x - lp };
		}

	Sample process_sample_and_mix( Sample x, Frequency cutoff_unwarped, Mix_1pole mix )
		{
		const auto filtered = process_sample( x, cutoff_unwarped );
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
	const auto cutoff_sampled = me.sample_function_over_domain( cutoff );

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
		[&]( Second t )
			{
			const Pole pole_c = pole( t );
			const Magnitude pole_abs = std::abs( pole_c );
			return std::make_pair( pole_abs, -pole_c.real() / pole_abs ); 
			},
		i );
	}

Audio base_filter_2pole_selector(
	const Audio & me,
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping,
	size_t i
	)
	{
	return base_filter_2pole_selector_single_func( me, [&]( Second t ){ return std::make_pair( cutoff( t ), damping( t ) ); }, i );
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

Audio Audio::filter_1pole_repeat(
	const Function<Second, Frequency> & cutoff,
	const uint16_t repeats
	) const
	{
	if( is_null() ) return Audio::create_null();
	
	Audio out( get_format() );

	auto cutoff_sampled = sample_function_over_domain( cutoff );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		std::vector<Filter_1Pole> filter_1poles( repeats, get_sample_rate() );
		for( Frame frame = 0; frame < get_num_frames(); ++frame )
			{
			for( Index i = 0; i < repeats; ++i )
				{
				const Audio & source = i == 0? *this : out;
				out.get_sample( channel, frame ) = filter_1poles[i].process_sample( source.get_sample( channel, frame ), cutoff_sampled[frame] )[0];
				}
			}
		}

	return out;
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

	const auto cutoff_sampled = me.sample_function_over_domain( cutoff );

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
	const auto cutoff_sampled = sample_function_over_domain( cutoff );

	std::vector<Audio> outs;
	outs.emplace_back( get_format() );
	outs.emplace_back( get_format() );
	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		Filter_1Pole filter( get_sample_rate() );
		for( Frame frame = 0; frame < get_num_frames(); ++frame )
			{
			auto filter_output = filter.process_sample( get_sample( channel, frame ), cutoff_sampled[frame] );
			outs[0].get_sample( channel, frame ) = filter_output[0];
			outs[1].get_sample( channel, frame ) = filter_output[1];
			}
		}
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

	const auto cutoff_sampled = me.sample_function_over_domain( cutoff );
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
	tilt.modify_volume_in_place( [&]( Second t ){ return decibel_to_amplitude( gain( t ) / 2 ); } );
	return tilt;
	}

Audio Audio::filter_1pole_highshelf(
	const Function<Second, Frequency> & cutoff,
	const Function<Second, Decibel> & gain,
	uint16_t order
	) const
	{
	if( is_null() ) return Audio::create_null();
	Audio tilt = base_filter_1pole_butterworth_tilt( *this, order, cutoff, [&]( Second t ){ return -gain( t ); } );
	tilt.modify_volume_in_place( [&]( Second t ){ return decibel_to_amplitude( gain( t ) / 2 ); } );
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

	const auto cutoff_sampled = me.sample_function_over_domain( cutoff );
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
	const std::vector<Decibel> gain_sampled = me.sample_function_over_domain( gain );
	std::vector<Amplitude> Ms( gain_sampled.size() );
	std::transform( gain_sampled.begin(), gain_sampled.end(), Ms.begin(), [&]( Decibel dB ){ return decibel_to_amplitude( dB / ( 2 * order ) ); } );
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
		[&]( Second t ){ return gain( t ) / 2.0f; }, 
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
		[&]( Second t ){ return -gain( t ); }, 
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
		[&]( Second t ){ return gain( t ) / 2.0f; }, 
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

Audio filter_1pole_allpass_n(
	const Audio & me,
	const Function<Second, Frequency> & cutoff,
	uint16_t order
	)
	{
	const auto cutoff_sampled = me.sample_function_over_domain( cutoff );

	Audio allpass = Audio::create_from_format( me.get_format() );

	for( Channel channel = 0; channel < me.get_num_channels(); ++channel )
		{
		std::vector<Filter_1Pole> filters( order, me.get_sample_rate() );
		for( Frame frame = 0; frame < me.get_num_frames(); ++frame )
			{
			Sample allpassed;
			const Frequency w = cutoff_sampled[frame];
			for( Index filter_i = 0; filter_i < filters.size(); ++filter_i )
				{
				const Sample input_sample = filter_i == 0? me.get_sample( channel, frame ) : allpassed;
				allpassed = filters[filter_i].process_sample_and_mix( input_sample, w, { 1.0f, -1.0f } );
				}
			allpass.get_sample( channel, frame ) = allpassed;
			}
		}

	return allpass;
	}

Audio Audio::filter_1pole_multinotch(
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & wet_dry,
	uint16_t order,
	bool invert
	) const
	{
	if( is_null() ) return Audio::create_null();
	
	// See section 11.1
	Audio allpass = filter_1pole_allpass_n( *this, cutoff, order );
	if( invert ) allpass.modify_volume_in_place( -1 );
	Function<Second, Amplitude> wet_dry_inv = [&]( Second t ){ return 1.0f - wet_dry( t ); };
	return Audio::mix_variable_gain( { this, &allpass }, {}, { &wet_dry, &wet_dry_inv } );
	}

Audio filter_2pole_allpass_n(
	const Audio & me,
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping,
	uint16_t order
	)
	{
	const auto cutoff_sampled = me.sample_function_over_domain( cutoff );
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
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping,
	const Function<Second, float> & wet_dry,
	uint16_t order,
	bool invert
	) const
	{
	if( is_null() ) return Audio::create_null();

	// See section 11.2
	Audio allpass = filter_2pole_allpass_n( *this, cutoff, damping, order );
	if( invert ) allpass.modify_volume_in_place( -1 );
	Function<Second, Amplitude> wet_dry_inv = [&]( Second t ){ return 1.0f - wet_dry( t ); };
	return Audio::mix_variable_gain( { this, &allpass }, {}, { &wet_dry, &wet_dry_inv } );
	}

Audio Audio::filter_comb(
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & wet_dry,
	const Function<Second, float> & feedback,
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

	const auto cutoff_sampled = sample_function_over_domain( cutoff );
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
