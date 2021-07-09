#include "flan/Audio.h"

#include <iostream>
#include <algorithm>
#include <complex>
#include <numeric>

#include "flan/Spectrum.h"
#include "flan/WindowFunctions.h"
#include "WDL/resample.h"

static const float pi = std::acos( -1.0f );
static const float pi2 = pi * 2.0f;

using namespace flan;

Audio::Audio() 
	: AudioBuffer() 
	{}

Audio::Audio( const AudioBuffer::Format & other ) 
	: AudioBuffer( other ) 
	{}

Audio::Audio( const std::string & filename ) 
	: AudioBuffer( filename ) 
	{}

Audio::Audio( const AudioBuffer & other )
	: AudioBuffer( other )
	{
	}

//========================================================
// Helper Functions
//========================================================

static bool doSampleRatesMatch( const std::vector<Audio> & ins )
	{
	//Check if all sample rates match the first file
	uint32_t sampleRate = ins[0].getSampleRate();
	for( auto & in : ins ) 
		if( in.getSampleRate() != sampleRate )
			{
			std::cout << "Mismatched sample rates" << std::endl;
			return false;
			}
	return true;
	}

static bool doChannelCountsMatch( const std::vector<Audio> & ins )
	{
	if( ins.size() == 0 ) return true;

	//Check if all channel counts match the first one
	uint32_t numchannels = ins[0].getNumChannels();
	for( auto & in : ins ) 
		if( in.getNumChannels() != numchannels )
			{
			std::cout << "Mismatched channel count" << std::endl;
			return false;
			}
	return true;
	}

static Channel getMaxNumChannels( const std::vector<Audio> & ins )
	{
	return std::max_element( ins.begin(), ins.end(), []( const Audio & a, const Audio & b )
		{ 
		return a.getNumChannels() < b.getNumChannels();
		} )->getNumChannels();
	}

static Frame getMaxNumFrames( const std::vector<Audio> & ins )
	{
	return std::max_element( ins.begin(), ins.end(), []( const Audio & a, const Audio & b )
		{ 
		return a.getNumFrames() < b.getNumFrames();
		} )->getNumFrames();
	}

//========================================================
// Overloads
//========================================================

Audio Audio::operator+( const Audio & other ) const
	{
	return mix( { *this, other } );
	}

Audio Audio::operator-() const	
	{
	return invertPhase();
	}

Audio Audio::operator-( const Audio & other ) const
	{
	return *this + -other;
	}

//========================================================
// Procs
//========================================================

Sample Audio::getSampleInterpolated( uint32_t channel, float t, Interpolator interp ) const
	{
	const float q = std::floor( t );
	const float r = t - q;
	const Sample s1 = getSample( channel, q     );
	const Sample s2 = getSample( channel, q + 1 );
	return s1 + ( s2 - s1 ) * interp( r );
	}

Audio Audio::invertPhase( flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	Audio out = deepCopy();

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		for( Frame frame = 0; frame < getNumFrames(); ++frame )
			{
			flan_CANCEL_POINT( Audio() );
			out.setSample( channel, frame, -getSample( channel, frame ) );
			}

	return out;
	}

Audio Audio::modifyVolume( Func1x1 volumeLevel, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	Audio out( getFormat() );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		for( Frame frame = 0; frame < getNumFrames(); ++frame )
			{
			flan_CANCEL_POINT( Audio() );
			const float calculatedSample = getSample( channel, frame ) * volumeLevel( frameToTime() * frame );
			out.setSample( channel, frame, calculatedSample );
			}

	return out;
	}

Audio Audio::setVolume( Func1x1 level, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	// Divide by getMaxSampleMagnitude to normalize, multiply by level to set
	const Sample maxMag = getMaxSampleMagnitude();
	if( maxMag == 0 ) return *this;
	return modifyVolume( [ &level, maxMag ]( float t ){ return level(t) / maxMag; }, canceller );
	}

Audio Audio::shift( Time shift, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	if( shift <= 0 ) return cut( -shift, getLength(), canceller );

	const Frame shiftFrames = shift * timeToFrame();

	auto format = getFormat();
	format.numFrames += shiftFrames;
	Audio out( format );
	out.clearBuffer();

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		for( Frame frame = 0; frame < getNumFrames(); ++frame )
			{
			flan_CANCEL_POINT( Audio() );
			out.setSample( channel, frame + shiftFrames, getSample( channel, frame ) );
			}
	}

Audio Audio::waveshape( Func1x1 shaper, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	Audio out( getFormat() );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		for( Frame frame = 0; frame < getNumFrames(); ++frame )
			{
			flan_CANCEL_POINT( Audio() );
			out.setSample( channel, frame, shaper( getSample( channel, frame ) ) );
			}

	return out;
	}

Audio Audio::pan( Func1x1 panAmount, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	//Stereo panning algorithm
	auto stereoPan = [this, &canceller]( Func1x1 panAmount )
		{
		Audio out( getFormat() );

		static const float sqrt2 = std::sqrt( 2.0f );

		// Channel 0
		for( Frame frame = 0; frame < getNumFrames(); ++frame )
			{
			flan_CANCEL_POINT( Audio() );
			const float pan = std::clamp( panAmount( frameToTime() * frame ), -1.0f, 1.0f );
			const float angle = pi / 4.0f * ( pan + 3.0f ); // pi/2 <= angle <= pi 
			const float scale = sin( angle ) * sqrt2;
			out.setSample( 0, frame, getSample( 0, frame ) * scale );
			}
		// Channel 1
		for( Frame frame = 0; frame < getNumFrames(); ++frame )
			{
			flan_CANCEL_POINT( Audio() );
			const float pan = std::clamp( panAmount( frameToTime() * frame ), -1.0f, 1.0f );
			const float angle = pi / 4.0f * ( pan + 1.0f ); // 0 <= angle <= pi/2
			const float scale = sin( angle ) * sqrt2;
			out.setSample( 1, frame, getSample( 1, frame ) * scale );
			}

		return out;
		};
	
	switch( getNumChannels() )
		{
		case 1: //Panning mono: cast to stereo and do stereo pan
			return convertToStereo().pan( panAmount );
			break;
		case 2: //Panning stereo: Use stereo pan as defined above
			return stereoPan( panAmount );
			break;
		default:
			std::cout << "I don't know how to pan that number of channels" << std::endl;
			return Audio();
		}
	}

Audio Audio::widen( Func1x1 widenAmount, flan_CANCEL_ARG_CPP ) const
	{
	flan_FUNCTION_LOG
	return convertToMidSide().pan( widenAmount, canceller ).convertToLeftRight( canceller );
	}

Audio Audio::reverse( flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	Audio out( getFormat() );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		for( Frame sample = 0; sample < getNumFrames(); ++sample )
			{
			flan_CANCEL_POINT( Audio() );
			out.setSample( channel, sample, getSample( channel, getNumFrames() - 1 - sample ) );
			}

	return out;
	}

Audio Audio::cut( Time startTime, Time endTime, Time startFade, Time endFade, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	const Frame startFrame = std::clamp( Frame( timeToFrame() * startTime ), 0, getNumFrames() - 1 );
	const Frame endFrame   = std::clamp( Frame( timeToFrame() * endTime   ), 0, getNumFrames() - 1 );
	const Frame startFadeFrames = startFade * timeToFrame();
	const Frame endFadeFrames = endFade * timeToFrame();
	
	return cutFrames( startFrame, endFrame, startFadeFrames, endFadeFrames );
	}

Audio Audio::cutFrames( Frame start, Frame end, Frame startFade, Frame endFade, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	// Input validation
	if( end <= start ) return Audio();
	// Any fade errors are validated in Audio::fades
	
	auto format = getFormat();
	format.numFrames = end - start;
	Audio out( format );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		for( Frame frame = 0; frame < out.getNumFrames(); ++frame )
			{
			flan_CANCEL_POINT( Audio() );
			out.setSample( channel, frame, getSample( channel, start + frame ) );
			}

	return out.fadeFrames( startFade, endFade, Interpolators::sqrt, canceller );
	}

Audio Audio::repitch( Func1x1 factor, Time granulTime, uint32_t qual, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	Frame granul = granulTime * timeToFrame();

	// Input validation
	auto safeFactorInv = [this, &factor]( float count )
		{ 
		constexpr float bound = 1000.0f;
		const float factorEval = factor( float(count) / getSampleRate() );
		return std::clamp( 1.0f / factorEval, 1.0f / bound, bound ); 
		};
	if( granul == 0 ) granul = 1;
	if( qual > 2 )
		{
		std::cout << "Invalid quality in Audio::repitch" << std::endl;
		return Audio();
		}

	// Estimate output length
	float acclen = 0.0;
	for( int count = 0; count < getNumFrames(); count += granul )
		acclen += granul * safeFactorInv( count );

	auto format = getFormat();
	format.numFrames = std::ceil( acclen );
	Audio out( format );
	
	WDL_Resampler rs;
		 if( qual == 0 ) rs.SetMode( true,  0, true, 64 ); // reasonable sinc resampling
	else if( qual == 1 ) rs.SetMode( true,  1, false    ); // linear interpolation with simple filter
	else if( qual == 2 ) rs.SetMode( false, 0, false    ); // no interpolation, useful for dirty sound
	else 
		{
		std::cout << "Unknown quality in Audio::repitch: " << qual << std::endl;
		return Audio();
		}
	const Channel chs = getNumChannels();
	const Frame inFrames = getNumFrames();
	std::vector<Sample> rsoutbuf( chs * granul );

	Frame inFrame = 0;
	Frame outFrame = 0;
	while( inFrame < getNumFrames() )
		{
		flan_CANCEL_POINT( Audio() );
		const double factor_to_use = safeFactorInv( inFrame );
		rs.SetRates( getSampleRate(), double( getSampleRate() ) * factor_to_use );
		WDL_ResampleSample* rsinbuf = nullptr;
		const Frame wanted = rs.ResamplePrepare( granul, chs, &rsinbuf );

		for( Channel channel = 0; channel < chs; ++channel )
			for( Frame j = 0; j < wanted; ++j )
				{
				if( inFrame + j < inFrames )
					rsinbuf[ j * chs + channel ] = getSample( channel, inFrame + j );
				else
					rsinbuf[ j * chs + channel ] = 0.0;
				}

		// Process
		rs.ResampleOut( rsoutbuf.data(), wanted, granul, chs );

		// Copy resampled data to output buffer
		for( Channel i = 0; i < chs; ++i )
			for( Frame j = 0; j < granul; ++j )
				if( outFrame + j < acclen )
					out.setSample( i, outFrame + j, rsoutbuf[ j  *chs + i ] );

		outFrame += granul;
		inFrame += wanted;
		}

	return out;
	}

Audio Audio::convolve( const std::vector<Func1x1> & ir, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	auto format = getFormat();
	format.numFrames = getNumFrames() + ir.size();
	Audio out( format );

	for( Channel channel = 0; channel < out.getNumChannels(); ++channel )
		for( Frame outFrame = 0; outFrame < out.getNumFrames(); ++outFrame )
			{
			Sample convolutionSum = 0;
			for( Frame irFrame = 0; irFrame < ir.size(); ++irFrame )
				{
				flan_CANCEL_POINT( Audio() );
				const Frame inFrame = outFrame + 1 - irFrame;
				if( 0 <= inFrame && inFrame < getNumFrames() )
					convolutionSum += ir[irFrame]( out.frameToTime() * outFrame ) * getSample( channel, inFrame );
				}
			out.setSample( channel, outFrame, convolutionSum );
			}

	return out;
	}

Audio Audio::convolve( const Audio & ir, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );
	if( ir.isNull() ) return Audio();

	auto format = getFormat();
	format.numFrames = getNumFrames() + ir.getNumFrames();
	Audio out( format );

	for( Channel channel = 0; channel < out.getNumChannels(); ++channel )
		for( Frame outFrame = 0; outFrame < out.getNumFrames(); ++outFrame )
			{
			Sample convolutionSum = 0;
			for( Frame irFrame = 0; irFrame < ir.getNumFrames(); ++irFrame )
				{
				flan_CANCEL_POINT( Audio() );
				const Frame inFrame = outFrame + 1 - irFrame;
				if( 0 <= inFrame && inFrame < getNumFrames() )
					convolutionSum += ir.getSample( channel, irFrame ) * getSample( channel, inFrame );
				}
			out.setSample( channel, outFrame, convolutionSum );
			}

	return out;
	}

Audio Audio::iterate( uint32_t n, Audio::Mod mod, bool feedback, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );
	// Half a length is removed to avoid length rounding changing the number of events generated
	return texture( getLength() * ( float( n ) - 0.5f ), 1.0f / getLength(), 0, mod, feedback, canceller );
	}

Audio Audio::delay( Time length, Func1x1 delayTime, Func1x1 decay, Audio::Mod mod, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	if( length <= 0 ) return Audio();
	return texture( length, 
		[&delayTime, this]( float t )
			{ 
			if( delayTime( t ) <= 0 ) return 1.0f / float( getSampleRate() );
			return 1.0f / delayTime( t );
			},
		0, [&decay, &mod, &canceller]( const Audio & in, float t )
			{ 
			if( t == 0 ) 
				return in;
			const float currentDecay = decay( t );
			return Audio( mod ? mod( in, t ) : in ).modifyVolume( currentDecay, canceller );
			},
		true, canceller );
	}

Audio Audio::fade( Time start, Time end, Interpolator interp, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	const Frame startFrames = start * timeToFrame();
	const Frame endFrames = end * timeToFrame();

	return fadeFrames( startFrames, endFrames, interp, canceller );
	}

Audio Audio::fadeFrames( Frame start, Frame end, Interpolator interp, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	// Input validation
	start = std::max( 0, start );
	end   = std::max( 0, end   );
	if( start + end > getNumFrames() )
		{
		const float scale = float( getNumFrames() ) / ( start + end );
		start = std::floor( start * scale );
		end	= std::floor( end * scale );
		}
	if( start == 0 && end == 0 ) return *this;

	Audio out = deepCopy();

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		for( Frame frame = 0; frame < start; ++frame )
			{
			flan_CANCEL_POINT( Audio() );
			const float fadeAmt = interp( float( frame ) / start );
			out.getSample( channel, frame ) *= fadeAmt;
			}
		for( Frame frame = 0; frame < end; ++frame )
			{
			flan_CANCEL_POINT( Audio() );
			const float fadeAmt = interp( float( frame ) / end );
			out.getSample( channel, getNumFrames() - 1 - frame ) *= fadeAmt;
			}
		}

	return out;
	}

Audio Audio::lowPass( Func1x1 cutoff, uint32_t taps, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	if( taps == 0 ) return *this;

	std::vector<Func1x1> ir;
	for( uint32_t tap = 0; tap <= taps; ++tap ) 
		ir.push_back( [&cutoff, this, taps, tap]( float t )
			{
			const float n = float(tap) - int(taps) / 2.0f;
			const float N = getSampleRate();
			const float K = 2.0f * cutoff(t);
			if( n == 0 ) return K / N;
			else return std::sin( pi * n * K / N ) / ( N * sin( pi * n / N ) ) 
				* Windows::Hann( float(tap) / float(taps) ); 
			} );

	return convolve( ir, canceller );
	}

Audio Audio::grainSelect( Time length, Func1x1 grainsPerSecond, Func1x1 scatter, 
		Func1x1 selection, Func1x1 grainLength, Func1x1 fade, Mod mod, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	// Input validation
	if( length <= 0 ) return Audio();

	return texture( length, grainsPerSecond, scatter, 
		[&selection, &grainLength, &fade, &mod, &canceller]( const Audio & in, float t )
			{
			const Time selection_c = selection( t );
			const Time grainLength_c = grainLength( t );
			const Time fade_c = fade( t );
			const Audio grain = in.cut( selection_c, selection_c + grainLength_c, fade_c, fade_c, canceller );
			return mod? mod( grain, t ) : grain;
			}, 
		false, canceller );
	}

std::vector<Audio> Audio::chop( Time sliceLength, Time fade, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( std::vector<Audio>() );

	// Input validation
	if( sliceLength <= 0 ) return std::vector<Audio>();

	// Get all but last slice
	std::vector<Audio> outs;
	Time slicePosition = 0;
	do 
		{
		outs.emplace_back( cut( slicePosition, slicePosition + sliceLength, fade, fade, canceller ) );
		slicePosition += sliceLength;
		}
		while( slicePosition < getLength() );

	return outs;
	}

Audio Audio::rearrange( Time sliceLength, Time fade, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( Audio() );

	// + fade here accounts for + fade / 2 at both ends
	auto chops = chop( sliceLength + fade, fade, canceller );
	if( chops.size() < 2 )
		return Audio();

	chops.pop_back(); // Final slice usually isn't the correct length

	std::random_device rd;
	std::mt19937 g( rd() );
	std::shuffle( chops.begin(), chops.end(), g );

	return Audio::join( chops, fade, canceller );
	}
	

//========================================================
// Multi-In Procs
//========================================================

Audio Audio::mix( const std::vector<Audio> & ins, std::vector<Func1x1> balances, 
	std::vector<Time> startTimes, flan_CANCEL_ARG_CPP )
	{
	flan_FUNCTION_LOG;

	// Input validation
	if( ins.size() == 0 ) return Audio();
	if( balances.size() < ins.size() ) 
		balances.insert( balances.end(), ins.size() - balances.size(), 1.0f );
	if( startTimes.size() < ins.size() )
		startTimes.insert( startTimes.end(), ins.size() - startTimes.size(), 0.0f );

	const float outTimeToFrame = ins[0].timeToFrame();

	// Output setup
	auto format = ins[0].getFormat();
	format.numChannels = getMaxNumChannels( ins );
	//Find how long the longest thing being mixed will last
	for( uint32_t in = 0; in < ins.size(); ++in )
		{
		format.numFrames = std::max( 
			format.numFrames, 
			ins[in].getNumFrames() + Frame( startTimes[in] * outTimeToFrame ) );
		}
	Audio out( format );
	out.clearBuffer();

	// For each input, write to output
	for( uint32_t in = 0; in < ins.size(); ++in )
		for( Channel channel = 0; channel < ins[in].getNumChannels(); ++channel )
			for( Frame frame = 0; frame < ins[in].getNumFrames(); ++frame )
				{
				flan_CANCEL_POINT( Audio() );
				const Frame outFrame = frame + startTimes[in] * outTimeToFrame;
				const Sample sample = ins[in].getSample( channel, frame ) 
					* balances[in]( ins[0].frameToTime() * outFrame );
				if( 0 <= outFrame && outFrame < format.numFrames )
					out.getSample( channel, outFrame ) += sample;
				}

	return out;
	}

Audio Audio::join( const std::vector<Audio> & ins, float fade, flan_CANCEL_ARG_CPP )
	{
	flan_FUNCTION_LOG;

	if( ins.empty() ) return Audio();
	fade = std::max( 0.0f, fade );

	// Get input Audio lengths
	std::vector<Time> jumps( ins.size() );
	jumps[0] = 0;
	std::transform( ins.begin(), ins.end() - 1, jumps.begin() + 1, []( const Audio & a ){ return a.getLength(); } );

	// Sum lengths to get mix positions
	std::vector<float> times( ins.size() );
	std::partial_sum( jumps.begin(), jumps.end(), times.begin(), [fade]( float a, float b ){ return a + b - fade / 2.0f; } );

	return mix( ins, std::vector<Func1x1>(), times, canceller );
	}

Audio Audio::select( const std::vector<Audio> & ins, Func1x1 selection, std::vector<Time> startTimes, flan_CANCEL_ARG_CPP )
	{
	flan_FUNCTION_LOG;

	// Generate balances from selection
	std::vector<Func1x1> balances;
	for( int i = 0; i < ins.size(); ++i )
		{
		balances.push_back( [selection, i]( Time t )
			{
			const float distance = std::abs( selection( t ) - i );
			if( distance >= 1 ) return 0.0f;
			else return std::sqrt( 1.0f - distance );
			} );
		}
		
	return mix( ins, balances, startTimes, canceller );
	}
