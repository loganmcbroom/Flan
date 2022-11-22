#include "flan/PVOC.h"

#include <iostream>
#include <complex>
#include <array>
#include <numeric>

namespace flan {

static const float pi = std::acos( -1.0f );

PVOC PVOC::getFrame( Time time ) const
	{
	flan_PROCESS_START( PVOC() );

	// getBinInterpolated does bounds checking
	const float selectedFrame = time * timeToFrame();

	auto format = getFormat();
	format.numFrames = 1;
	PVOC out( format );

	for( Channel channel = 0; channel < out.getNumChannels(); ++channel )
		for( Bin bin = 0; bin < out.getNumBins(); ++bin )
			out.setMF( channel, 0, bin, getBinInterpolated( channel, selectedFrame, bin ) );

	return out;
	}

MF PVOC::getBinInterpolated( Channel channel, float frame, float bin, Interpolator i ) const
	{
	frame = std::clamp( frame, 0.0f, float( getNumFrames() - 1 ) );
	bin   = std::clamp( bin,   0.0f, float( getNumBins()   - 1 ) );

	const std::array< MF, 4 > p = 
		{
		getMF( channel, floor( frame ), floor( bin ) ),
		getMF( channel, ceil ( frame ), floor( bin ) ),
		getMF( channel, ceil ( frame ), ceil ( bin ) ),
		getMF( channel, floor( frame ), ceil ( bin ) )
		};

	const float l = i( frame - floor( frame ) );
	const float m = i( bin   - floor( bin   ) );
	const float mI = 1.0f - m;
	const float lI = 1.0f - l;

	return {
		   mI * ( lI * p[0].m + l * p[1].m ) + m * ( lI * p[3].m + l * p[2].m ),
		   mI * ( lI * p[0].f + l * p[1].f ) + m * ( lI * p[3].f + l * p[2].f )
		   };
	}

MF PVOC::getBinInterpolated( Channel channel, float frame, Bin bin, Interpolator i ) const
	{
	frame = std::clamp( frame, 0.0f, float( getNumFrames() - 1 ) );

	const MF & l = getMF( channel, floor( frame ), bin );
	const MF & h = getMF( channel, ceil ( frame ), bin );

	const float mix = i( frame - floor( frame ) );

	return { 
		   ( 1.0f - mix ) * l.m + mix * h.m,
		   ( 1.0f - mix ) * l.f + mix * h.f
		   };
	}

MF PVOC::getBinInterpolated( Channel channel, Frame frame, float bin, Interpolator i ) const
	{
	bin = std::clamp( bin, 0.0f, float( getNumBins() - 1 ) );

	const auto l = getMF( channel, frame, floor( bin ) );
	const auto h = getMF( channel, frame, ceil ( bin ) );

	const float mix = i( bin - floor( bin ) );

	return  { 
			( 1.0f - mix ) * l.m + mix * h.m,
			( 1.0f - mix ) * l.f + mix * h.f
			};
	}

//========================================================================
//	Selection
//========================================================================

PVOC PVOC::select( Time length, Func2x2 selector, Interpolator interp, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( PVOC() );

	// Input validation
	if( length <= 0 ) return PVOC();

	auto format = getFormat();
	format.numFrames = timeToFrame() * length;
	PVOC out( format );

	for( Channel channel = 0; channel < out.getNumChannels(); ++channel )
		for( Frame frame = 0; frame < out.getNumFrames(); ++frame )
			{
			for( Bin bin = 0; bin < getNumBins(); ++bin )
				{
				flan_CANCEL_POINT( PVOC() );
				//No need to clamp, getBinInterpolated does bounds checking
				vec2 selectedPoint = selector( frame * frameToTime(), bin * binToFrequency() );
				selectedPoint.x() *= timeToFrame();
				selectedPoint.y() *= frequencyToBin();
				out.setMF( channel, frame, bin,  
					getBinInterpolated( channel, selectedPoint.x(), selectedPoint.y(), interp ) );
				}
			}

	return out;
	}

PVOC PVOC::freeze( const std::vector<std::array<Time,2>> & timing, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( PVOC() );

	using TimingPair = std::array<Frame,2>;

	//Convert times into frames
	std::vector<TimingPair> timingFrames( timing.size() );
	std::transform( timing.begin(), timing.end(), timingFrames.begin(), [this]( const std::array<Time,2> & i )
		{ 
		return TimingPair{ 
			std::clamp( Frame( i[0] * timeToFrame() ), 0, Frame( getNumFrames() - 1 ) ), 
			std::max( Frame( i[1] * timeToFrame() ), 0 ) };
		});

	//Sort by occurance frame
	std::sort( timingFrames.begin(), timingFrames.end(), []( TimingPair & a, TimingPair & b ){ return a[0] < b[0]; } );

	// Remove simultaneous events
	auto timingEnd = std::unique( timingFrames.begin(), timingFrames.end(), 
		[]( const TimingPair & a, const TimingPair & b )
			{
			return a[0] == b[0];
			} );
	timingFrames.erase( timingEnd, timingFrames.end() ); 

	float totalFreezeFrames = 0;
	for( auto & i : timingFrames ) totalFreezeFrames += i[1];

	auto format = getFormat();
	format.numFrames += totalFreezeFrames;
	PVOC out( format );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		uint32_t timingIndex = 0;
		for( Frame inFrame = 0, outFrame = 0; inFrame < getNumFrames(); ++inFrame )
			{
			//If it's time to freeze, freeze
			if( inFrame == timingFrames[timingIndex][0] )
				{
				for( Frame freezeFrame = 0; freezeFrame < timingFrames[timingIndex][1]; ++freezeFrame )
					{
					flan_CANCEL_POINT( PVOC() );
					for( Bin bin = 0; bin < getNumBins(); ++bin )
						out.setMF( channel, outFrame, bin, getMF( channel, inFrame, bin ) );
					++outFrame;
					}
				}
			else //It's not time to freeze
				{
				flan_CANCEL_POINT( PVOC() );
				for( Bin bin = 0; bin < getNumBins(); ++bin )
					out.setMF( channel, outFrame, bin, getMF( channel, inFrame, bin ) );
				++outFrame;
				}
			}
		}

	return out;
	}


//========================================================================
// Combinations
//========================================================================

PVOC PVOC::replaceAmplitudes( const PVOC & ampSource, Func2x1 amount, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( PVOC() );

	// Input validation
	if( ampSource.isNull() ) return PVOC();
	Func2x1 amount_s = Func2x1::clamp( amount, 0.0f, 1.0f );

	PVOC out( getFormat() );
	out.clearBuffer();

	const uint32_t numChannels = std::min( ampSource.getNumChannels(), getNumChannels() );
	const uint32_t numFrames   = std::min( ampSource.getNumFrames(),   getNumFrames()   );
	const uint32_t numBins     = std::min( ampSource.getNumBins(),     getNumBins()     );

	for( Channel channel = 0; channel < numChannels; ++channel )
		for( Frame frame = 0; frame < numFrames; ++frame )
			for( Bin bin = 0; bin < numBins; ++bin )
				{
				flan_CANCEL_POINT( PVOC() );
				const auto currentBin = getMF( channel, frame, bin );
				const float amount_c = amount_s( frame * frameToTime(), bin * binToFrequency() );
				out.setMF( channel, frame, bin,
					{
					float( ampSource.getMF( channel, frame, bin ).m * amount_c + currentBin.m * ( 1.0f - amount_c ) ),
					float( currentBin.f )
					} );
				}

	return out;
	}

PVOC PVOC::subtractAmplitudes( const PVOC & ampSource, Func2x1 amount, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( PVOC() );

	// Input validation
	if( ampSource.isNull() ) return PVOC();

	PVOC out = deepCopy();

	const uint32_t numChannels = std::min( ampSource.getNumChannels(), getNumChannels() );
	const uint32_t numFrames   = std::min( ampSource.getNumFrames(),   getNumFrames()   );
	const uint32_t numBins     = std::min( ampSource.getNumBins(),     getNumBins()     );

	for( Channel channel = 0; channel < numChannels; ++channel )
		for( Frame frame = 0; frame < numFrames; ++frame )
			for( Bin bin = 0; bin < numBins; ++bin )
				{
				flan_CANCEL_POINT( PVOC() );
				auto & outBin = out.getMF( channel, frame, bin );
				const float amount_c = amount( frame * frameToTime(), bin * binToFrequency() );
				outBin.m = std::abs( outBin.m - ampSource.getMF( channel, frame, bin ).m * amount_c );
				}

	return out;
	}


//========================================================================
// Generation
//========================================================================

PVOC PVOC::synthesize( Time length, Func1x1 freq, Func2x1 harmonicWeights, flan_CANCEL_ARG_CPP )
	{
	PVOC::Format format;
	format.hopSize = 128;
	format.numBins = 2049;
	format.numChannels = 1;
	format.numFrames = length * 48000 / 128; // Sample rate over hop size
	format.sampleRate = 48000;
	format.windowSize = 2048;

	PVOC out( format );
	out.clearBuffer();
	const float scale = std::sqrt( out.getDFTSize() ); // Don't think too much about this constant

	for( Frame frame = 0; frame < out.getNumFrames(); ++frame )
		{
		const Frequency freq_c = freq( frame * out.frameToTime() );
		if( freq_c <= 1.0 ) continue;
		for( int harmonic = 0; harmonic < std::floor( out.getHeight() / freq_c ); ++harmonic )
			{
			const float mag = harmonicWeights( frame * out.frameToTime(), harmonic ) * scale;
			const Frequency harmonicFreq = freq_c * ( harmonic + 1 );
			const Bin bin = harmonicFreq * out.frequencyToBin();
			out.getMF( 0, frame, bin ) = { mag, harmonicFreq };
			}
		}
	
	return out;
	}

//========================================================================
// Uncategorized
//========================================================================

PVOC harmonicScaler( const PVOC & me, Func2x1 series, std::function< Frequency ( Frequency, int )> harmonicFunc, int numSamples, flan_CANCEL_ARG_CPP )
	{
	PVOC out( me.getFormat() );
	out.clearBuffer();

	for( Channel channel = 0; channel < me.getNumChannels(); ++channel )
		for( Frame frame = 0; frame < me.getNumFrames(); ++frame )
			{
			// Samples for this frame - note this could take the time-freq pos and an index for added flexibility. I don't think it's needed.
			std::vector<float> seriesSamples( numSamples );
			for( int i = 0; i < seriesSamples.size(); ++i )
				seriesSamples[i] = series( frame * me.frameToTime(), i );

			for( Bin bin = 0; bin < me.getNumBins(); ++bin )
				{
				flan_CANCEL_POINT( PVOC() );
				const MF & source = me.getMF( channel, frame, bin );
				if( source.f <= 1.0f ) continue;
				int i = 0;
				for( Frequency freq = source.f; freq < me.getHeight() && i < seriesSamples.size(); freq = harmonicFunc( source.f, i ) )
					{
					MF & dest = out.getMF( channel, frame, freq * me.frequencyToBin() );
					if( dest.m < source.m * seriesSamples[i] )
						{
						dest.m = source.m * seriesSamples[i];
						dest.f = freq;
						}
					++i;
					}
				}
			}

	return out;
	}

PVOC PVOC::addOctaves( Func2x1 series, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( PVOC() );
	return harmonicScaler( *this, series, []( Frequency f, int i ){ return f * std::pow( 2, i ); }, std::ceil( std::log2( getHeight() ) ), canceller );
	}

PVOC PVOC::addHarmonics( Func2x1 series, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( PVOC() );
	return harmonicScaler( *this, series, []( Frequency f, int i ){ return f * i; }, getNumBins(), canceller );
	}

PVOC PVOC::shape( Func2x2 shaper, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( PVOC() );

	auto shaperConverter = [&shaper]( const MF & mf )
		{ 
		const vec2 out = shaper( mf.m, mf.f );
		return MF{ out.x(), out.y() };
		};

	PVOC out( getFormat() );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		for( Frame frame = 0; frame < getNumFrames(); ++frame )
			for( Bin bin = 0; bin < getNumBins(); ++bin )
				{
				flan_CANCEL_POINT( PVOC() );
				const MF & currentBin = getMF( channel, frame, bin );
				out.setMF( channel, frame, bin, shaperConverter( currentBin ) );
				}
			
	return out;
	}

PVOC PVOC::perturb( Func2x1 magSigma, Func2x1 frqSigma, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( PVOC() );

	// Input validation
	Func2x1 magSigma_s = Func2x1::max( magSigma, 0 );
	Func2x1 frqSigma_s = Func2x1::max( frqSigma, 0 );

	//Func2x1 magLogScalar = [this]( float, float f ){ return log2( 1.0f + f * frequencyToBin() ); };
	Func2x1 magPert = Func2x1::normalDistribution( 0, magSigma_s );
	Func2x1 frqPert = Func2x1::normalDistribution( 0, frqSigma_s );

	PVOC out( getFormat() );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		for( Frame frame = 0; frame < getNumFrames(); ++frame )
			{
			const float time = frame * frameToTime();
			for( Bin bin = 0; bin < getNumBins(); ++bin )
				{
				flan_CANCEL_POINT( PVOC() );
				const float freq = bin * binToFrequency();
				const MF & currentBin = getMF( channel, frame, bin );
				out.setMF( channel, frame, bin,
					{
					currentBin.m + magPert( time, freq ),
					currentBin.f + frqPert( time, freq )
					} );
				}
			}

	return out;
	}

PVOC predicateNLoudestPartials( const PVOC & me, Func1x1 numBins, 
	std::function< bool (uint32_t, uint32_t) > predicate, flan_CANCEL_ARG_CPP )
	{
	flan_FUNCTION_LOG;

	// Input validation
	Func1x1 numBins_s = Func1x1::clamp( numBins, 0, int( me.getNumBins() ) );

	PVOC out( me.getFormat() );
	out.clearBuffer();

	typedef std::pair<uint32_t, float> indexVolume;
	std::vector<indexVolume> indexAndVolumes( me.getNumBins() );

	for( Channel channel = 0; channel < me.getNumChannels(); ++channel )
		for( Frame frame = 0; frame < me.getNumFrames(); ++frame )
			{
			for( Bin bin = 0; bin < me.getNumBins(); ++bin )
				{
				flan_CANCEL_POINT( PVOC() );
				indexAndVolumes[bin].first = bin;
				indexAndVolumes[bin].second = me.getMF( channel, frame, bin ).m;
				}
			std::sort( indexAndVolumes.begin(), indexAndVolumes.end(), []( indexVolume& a, indexVolume& b ){ return abs(a.second) > abs(b.second); } );
			for( Bin bin = 0; bin < me.getNumBins(); ++bin )
				{
				flan_CANCEL_POINT( PVOC() );
				const uint32_t actualBin = indexAndVolumes[bin].first;
				if( predicate( bin, numBins_s( frame ) ) )
					out.setMF( channel, frame, actualBin, me.getMF( channel, frame, actualBin ) );
				else
					out.setMF( channel, frame, actualBin, { 0.0, me.getMF( channel, frame, actualBin ).f } );
				}
			}

	return out;
	}

PVOC PVOC::retainNLoudestPartials( Func1x1 numBins, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( PVOC() );
	return predicateNLoudestPartials( *this, numBins, []( uint32_t a, uint32_t b ){ return a < b; }, canceller );
	}

PVOC PVOC::removeNLoudestPartials( Func1x1 numBins, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( PVOC() );
	return predicateNLoudestPartials( *this, numBins, []( uint32_t a, uint32_t b ){ return a >= b; }, canceller );
	}

PVOC PVOC::resonate( Time length, Func2x1 decay, flan_CANCEL_ARG_CPP ) const
	{
	flan_PROCESS_START( PVOC() );

	// Input validation
	if( length < 0 ) 
		length = 0;
	Func2x1 decay_s = Func2x1::max( 0, decay );

	auto format = getFormat();
	format.numFrames = getNumFrames() + ceil( timeToFrame() * length );
	PVOC out( format );

	//Copy first frame into out
	for( Channel channel = 0; channel < out.getNumChannels(); ++channel )
		for( Bin bin = 0; bin < out.getNumBins(); ++bin )
			{
			flan_CANCEL_POINT( PVOC() );
			out.setMF( channel, 0, bin, getMF( channel, 0, bin ) );
			}
	
	const float frameToTimeConst = frameToTime();
	const float secondsPerFrame = frameToTime();
	const float binToFreqConst = binToFrequency();

	for( Channel channel = 0; channel < out.getNumChannels(); ++channel )
		for( Frame frame = 1; frame < out.getNumFrames(); ++frame )
			for( Bin bin = 0; bin < out.getNumBins(); ++bin )
				{
				flan_CANCEL_POINT( PVOC() );
				const float currentDecay = std::pow( decay_s( frameToTimeConst * frame, binToFreqConst * bin  ), secondsPerFrame );
				const float decayedAmp = out.getMF( channel, frame - 1, bin ).m * currentDecay;
				if( frame < getNumFrames() && getMF( channel, frame, bin ).m > decayedAmp )
					out.setMF( channel, frame, bin, getMF( channel, frame, bin ) );
				else
					out.setMF( channel, frame, bin, { decayedAmp, out.getMF( channel, frame - 1, bin ).f } );
				}

	return out;
	}

} //End namespace flan
