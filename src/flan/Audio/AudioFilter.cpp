#include "flan/Audio/Audio.h"

/*
https://ia601900.us.archive.org/5/items/the-art-of-va-filter-design-rev.-2.1.2/VAFilterDesign_2.1.2.pdf#chapter.10
*/

using namespace flan;

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

std::vector<Audio> Audio::filter_1pole_multimode( 
	const Function<Second, Frequency> & cutoff,
	bool return_low,
	bool return_high
	) const
	{

	/*
	See section 3.10 for the implementation.
	*/

	Audio low( return_low ? getFormat() : Audio::Format() );
	Audio high( return_high ? getFormat() : Audio::Format() );

	// Note the factor of 2pi. The reference book uses a non-standard Forier transform definition, much to my annoyance. This fixes that.
	// Saving a few cycles, the factor of 1/2 is baked into T, taking the 2pi down to pi.
	const Radian T_half = pi / getSampleRate(); 

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		float s = 0;
		for( Frame frame = 0; frame < getNumFrames(); ++frame )
			{
			const Sample x = getSample( channel, frame );

			const float w = prewarp( cutoff( frameToTime( frame ) ), T_half );
			
			const float g = w * T_half; 
			const float G = g / ( 1 + g );
			const float v = G * ( x - s );
			const float y_LP = v + s;
			s = y_LP + v;

			if( return_low ) low.getSample( channel, frame ) = y_LP;
			if( return_high ) high.getSample( channel, frame ) = x - y_LP;
			}
		}

	std::vector<Audio> outs;
	if( return_low )  outs.emplace_back( std::move( low ) );
	if( return_high ) outs.emplace_back( std::move( high ) );

	return outs;
	}

Audio Audio::filter_1pole_lowpass( 
	const Function<Second, Frequency> & cutoff 
	) const
	{
	return std::move( filter_1pole_multimode( cutoff, true, false )[0] );
	}

Audio Audio::filter_1pole_highpass( 
	const Function<Second, Frequency> & cutoff
	) const
	{
	return std::move( filter_1pole_multimode( cutoff, false, true )[0] );
	}

Audio Audio::filter_1pole_lowshelf( 
	const Function<Second, Frequency> & cutoff,
	const Function<Second, Decibel> & gain
	) const
	{

	/*
	See section 10.2 for the background for the shelving filters.
	For a 1-pole transfer function G_w(s)=w/(s+w), simplifying G(s/M)/G(MS) gives G_Mw(s) + M^2 * ( s / (s + w) ).
	That is, a cutoff scaled copy of a lowpass and an amplitude scaled copy of the corresponding highpass.
	This produces a high shelving filter which can then be scaled to obtain a low shelf.
	*/

	const auto M = [&]( Second t ){ return decibelToAmplitude( -gain( t ) / 2.0f ); };
	Audio lp = filter_1pole_lowpass( [&]( Second t ){ return M( t ) * cutoff( t ); } );
	const Audio hp = filter_1pole_highpass( cutoff );
	lp.modifyVolumeInPlace( [&]( Second t ){ return std::pow( M( t ), -2.0f ); } );

	return Audio::mix( std::vector<const Audio *>{ &lp, &hp } );
	}

Audio Audio::filter_1pole_highshelf( 
	const Function<Second, Frequency> & cutoff,
	const Function<Second, Decibel> & gain
	) const
	{
	const auto M = [&]( Second t ){ return decibelToAmplitude( gain( t ) / 2.0f ); };
	const Audio lp = filter_1pole_lowpass( [&]( Second t ){ return M( t ) * cutoff( t ); } );
	Audio hp = filter_1pole_highpass( cutoff );
	hp.modifyVolumeInPlace( [&]( Second t ){ return std::pow( M( t ), 2.0f ); } );

	return Audio::mix( std::vector<const Audio *>{ &lp, &hp } );
	}

std::vector<Audio> Audio::filter_svf_multimode(
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping,
	bool return_low,
	bool return_band,
	bool return_high
	) const
	{
	
	/*
	See section 4.4 for the implementation.
	*/

	Audio low ( return_low  ? getFormat() : Audio::Format() );
	Audio band( return_band ? getFormat() : Audio::Format() );
	Audio high( return_high ? getFormat() : Audio::Format() );

	const Radian T_half = pi / getSampleRate(); 

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		float s1 = 0;
		float s2 = 0;
		for( Frame frame = 0; frame < getNumFrames(); ++frame )
			{
			const Sample x = getSample( channel, frame );

			const Frequency w = prewarp( cutoff( frameToTime( frame ) ), T_half );
			const float R = damping( frameToTime( frame ) );

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

			if( return_low  ) low .getSample( channel, frame ) = lp;
			if( return_band ) band.getSample( channel, frame ) = bp * 2*R;
			if( return_high ) high.getSample( channel, frame ) = hp;
			}
		}

	std::vector<Audio> outs;
	if( return_low  ) outs.emplace_back( std::move( low  ) );
	if( return_band ) outs.emplace_back( std::move( band ) );
	if( return_high ) outs.emplace_back( std::move( high ) );

	return outs;
	}

Audio Audio::filter_svf_lowpass(
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping
	) const
	{
	return std::move( filter_svf_multimode( cutoff, damping, true, false, false )[0] );
	}

Audio Audio::filter_svf_bandpass(
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping
	) const
	{
	return std::move( filter_svf_multimode( cutoff, damping, false, true, false )[0] );
	}

Audio Audio::filter_svf_highpass(
	const Function<Second, Frequency> & cutoff,
	const Function<Second, float> & damping
	) const
	{
	return std::move( filter_svf_multimode( cutoff, damping, false, false, true )[0] );
	}

Audio filter_butterworth(
	const Audio & me,
	uint16_t order,
	const Function<Second, Frequency> & cutoff,
	bool lowpass
	)
	{
	/*
	See section 8.6 for details.
	*/

	auto svf_for_angle = [&]( const Audio & a, Radian theta )
		{ 
		const float R = std::cos( theta );
		return lowpass ? a.filter_svf_lowpass( cutoff, R ) : a.filter_svf_highpass( cutoff, R );
		};

	if( order == 0 ) return Audio();

	// The even/odd versions of this algorithm are similar enough to be factored further, but the code is easier to follow this way
	if( order % 2 == 0 )
		{
		const Radian theta = pi / ( 2.0f * order ); // Half the angular distance between poles for this order
		Audio out = svf_for_angle( me, theta );
		for( int i = 1; i < order / 2; ++i )
			out = svf_for_angle( out, theta * ( 1 + 2 * i ) );
		return out;
		}
	else
		{
		const Radian theta = pi / order; // Angular distance between poles for this order
		Audio out = lowpass ? me.filter_1pole_lowpass( cutoff ) : me.filter_1pole_highpass( cutoff );
		for( int i = 1; i < ( order + 1 ) / 2.0f; ++i )
			out = svf_for_angle( out, theta * i );
		return out;
		}	
	}

Audio Audio::filter_butterworth_lowpass(
	uint16_t order,
	const Function<Second, Frequency> & cutoff
	) const
	{
	return filter_butterworth( *this, order, cutoff, true );
	}

Audio Audio::filter_butterworth_highpass(
	uint16_t order,
	const Function<Second, Frequency> & cutoff
	) const
	{
	return filter_butterworth( *this, order, cutoff, false );
	}

Audio Audio::filter_butterworth_lowshelf(
	uint16_t order,
	const Function<Second, Frequency> & cutoff,
	const Function<Second, Decibel> & gain
	) const
	{
	const auto M = [&]( Second t ){ return decibelToAmplitude( -gain( t ) / 2.0f ); };
	Audio lp = filter_butterworth_lowpass( order, [&]( Second t ){ return M( t ) * cutoff( t ); } );
	const Audio hp = filter_butterworth_highpass( order, cutoff );
	lp.modifyVolumeInPlace( [&]( Second t ){ return std::pow( M( t ), -2.0f * order ); } );

	return Audio::mix( std::vector<const Audio *>{ &lp, &hp } );
	}