#include "flan/PV/PV.h"

#include <iostream>
#include <complex>
#include <array>
#include <numeric>
#include <ranges>

#include "flan/Utility/buffer_access.h"

using namespace std::ranges;

namespace flan {

PV PV::getFrame( Time time ) const
	{
	flan_PROCESS_START( PV() );

	const float selectedFrame = std::clamp( time * timeToFrame(), 0.0f, float( getNumFrames() - 1 ) );

	auto format = getFormat();
	format.numFrames = 1;
	PV out( format );

	for( Channel channel = 0; channel < out.getNumChannels(); ++channel )
		for( Bin bin = 0; bin < out.getNumBins(); ++bin )
			out.getMF( channel, 0, bin ) = getBinInterpolated( channel, selectedFrame, bin );

	return out;
	}

MF PV::getBinInterpolated( Channel channel, float frame, float bin, Interpolator i ) const
	{
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

MF PV::getBinInterpolated( Channel channel, float frame, Bin bin, Interpolator i ) const
	{
	const MF & l = getMF( channel, floor( frame ), bin );
	const MF & h = getMF( channel, ceil ( frame ), bin );

	const float mix = i( frame - floor( frame ) );

	return { 
		   ( 1.0f - mix ) * l.m + mix * h.m,
		   ( 1.0f - mix ) * l.f + mix * h.f
		   };
	}

MF PV::getBinInterpolated( Channel channel, Frame frame, float bin, Interpolator i ) const
	{
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

PV PV::select( Time length, const Func2x2 & selector, bool interpolateFrames, Interpolator interp ) const
	{
	flan_PROCESS_START( PV() );

	// Input validation
	if( length <= 0 ) return PV();

	auto format = getFormat();
	format.numFrames = timeToFrame() * length;
	PV out( format );

	auto selectorSamples = selector.sample( 0, out.getNumFrames(), frameToTime(), 0, getNumBins(), binToFrequency() );
	std::for_each( std::execution::par_unseq, selectorSamples.begin(), selectorSamples.end(), [&]( vec2 & v ){
		v.x() *= timeToFrame();
		v.y() *= frequencyToBin();
	 	} );

	for( Channel channel = 0; channel < out.getNumChannels(); ++channel )
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.getNumFrames() ), [&]( Frame frame )
			{
			for( Bin bin = 0; bin < out.getNumBins(); ++bin )
				{
				const vec2 s = selectorSamples[ buffer_access( bin, frame, getNumBins() ) ];

				if( s.x() < 0 || getNumFrames() - 1 <= s.x() ||
				    s.y() < 0 || getNumBins()   - 1 <= s.y() )
					continue;

				// Get the MFs surrounding the selected point s
				const vec2 sFloor = { std::floor( s.x() ), std::floor( s.y() ) };
				const std::array<MF, 4> MFs = {
					getMF( channel, sFloor.x()	  , sFloor.y()     ),
					getMF( channel, sFloor.x()	  , sFloor.y() + 1 ),
					getMF( channel, sFloor.x() + 1, sFloor.y()     ),
					getMF( channel, sFloor.x() + 1, sFloor.y() + 1 ),
					};
				
				// Get the weighted magnitudes of those points
				const vec2 mix = s - sFloor;
				const std::array<Magnitude, 4> w = {
					( 1.0f - mix.x() ) * ( 1.0f - mix.y() ) * MFs[0].m,
					(        mix.x() ) * ( 1.0f - mix.y() ) * MFs[1].m,
					(        mix.x() ) * (        mix.y() ) * MFs[2].m,
					( 1.0f - mix.x() ) * (        mix.y() ) * MFs[3].m };

				const auto maxWIter = std::max_element( w.begin(), w.end() );
				const int maxWIndex = std::distance( w.begin(), maxWIter );
				const MF selectedMF = MFs[maxWIndex];
				
				out.getMF( channel, frame, bin ) = selectedMF;
				}
			} );

	return out;
	}

PV PV::freeze( const std::vector<std::array<Time,2>> & timing ) const
	{
	flan_PROCESS_START( PV() );

	using TimingPair = std::array<Frame,2>;

	//Convert times into frames
	std::vector<TimingPair> timingFrames( timing.size() );
	std::ranges::transform( timing, timingFrames.begin(), [this]( const std::array<Time,2> & i )
		{ 
		return TimingPair{ 
			std::clamp( Frame( i[0] * timeToFrame() ), 0, Frame( getNumFrames() - 1 ) ), 
			std::max  ( Frame( i[1] * timeToFrame() ), 0 ) };
		});

	//Sort by occurance frame
	std::sort( std::execution::par_unseq, timingFrames.begin(), timingFrames.end(), []( TimingPair & a, TimingPair & b ){ return a[0] < b[0]; } );

	// Remove simultaneous events
	auto timingEnd = std::unique( timingFrames.begin(), timingFrames.end(), []( const TimingPair & a, const TimingPair & b ) {
		return a[0] == b[0];
		} );
	timingFrames.erase( timingEnd, timingFrames.end() ); 

	float totalFreezeFrames = 0;
	for( auto & i : timingFrames ) totalFreezeFrames += i[1];

	auto format = getFormat();
	format.numFrames += totalFreezeFrames;
	PV out( format );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		uint32_t timingIndex = 0;
		for( Frame inFrame = 0, outFrame = 0; inFrame < getNumFrames(); ++inFrame )
			{
			// If it's time to freeze, freeze
			if( inFrame == timingFrames[timingIndex][0] )
				{
				for( Frame freezeFrame = 0; freezeFrame < timingFrames[timingIndex][1]; ++freezeFrame )
					{
					for( Bin bin = 0; bin < getNumBins(); ++bin )
						out.setMF( channel, outFrame, bin, getMF( channel, inFrame, bin ) );
					++outFrame;
					}
				}
			else // It's not time to freeze
				{
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

PV PV::replaceAmplitudes( const PV & ampSource, const Func2x1 & amount ) const
	{
	flan_PROCESS_START( PV() );

	// Input validation
	if( ampSource.isNull() ) return PV();
	auto amountSamples = amount.sample( 0, getNumFrames(), frameToTime(), 0, getNumBins(), binToFrequency() );
	std::for_each( amountSamples.begin(), amountSamples.end(), []( float & amount ){ amount = std::clamp( amount, 0.0f, 1.0f ); } );

	PV out( getFormat() );
	out.clearBuffer();

	const uint32_t numChannels = std::min( ampSource.getNumChannels(), getNumChannels() );
	const uint32_t numFrames   = std::min( ampSource.getNumFrames(),   getNumFrames()   );
	const uint32_t numBins     = std::min( ampSource.getNumBins(),     getNumBins()     );

	for( Channel channel = 0; channel < numChannels; ++channel )
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( numFrames ), [&]( Frame frame )
			{
			for( Bin bin = 0; bin < numBins; ++bin )
				{
				const auto currentBin = getMF( channel, frame, bin );
				const float amount_c = amountSamples[buffer_access( bin, frame, getNumBins() )];
				out.setMF( channel, frame, bin,
					{
					float( ampSource.getMF( channel, frame, bin ).m * amount_c + currentBin.m * ( 1.0f - amount_c ) ),
					float( currentBin.f )
					} );
				}
			} );
	return out;
	}

PV PV::subtractAmplitudes( const PV & ampSource, const Func2x1 & amount ) const
	{
	flan_PROCESS_START( PV() );

	// Input validation
	if( ampSource.isNull() ) return PV();
	auto amountSamples = amount.sample( 0, getNumFrames(), frameToTime(), 0, getNumBins(), binToFrequency() );

	PV out = copy();

	const uint32_t numChannels = std::min( ampSource.getNumChannels(), getNumChannels() );
	const uint32_t numFrames   = std::min( ampSource.getNumFrames(),   getNumFrames()   );
	const uint32_t numBins     = std::min( ampSource.getNumBins(),     getNumBins()     );

	for( Channel channel = 0; channel < numChannels; ++channel )
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( numFrames ), [&]( Frame frame )
			{
			for( Bin bin = 0; bin < numBins; ++bin )
				{
				auto & outBin = out.getMF( channel, frame, bin );
				const float amount_c = amountSamples[buffer_access( bin, frame, getNumBins() )];
				outBin.m = std::abs( outBin.m - ampSource.getMF( channel, frame, bin ).m * amount_c );
				} 
			} );

	return out;
	}


//========================================================================
// Generation
//========================================================================

PV PV::synthesize( Time length, Func1x1 freq, const Func2x1 & harmonicWeights )
	{
	PV::Format format;
	format.hopSize = 128;
	format.numBins = 2049;
	format.numChannels = 1;
	format.sampleRate = 48000;
	format.numFrames = length * format.sampleRate / format.hopSize;
	format.windowSize = 2048;

	PV out( format );
	out.clearBuffer();
	const float scale = std::sqrt( out.getDFTSize() ); // Don't think too much about this constant

	auto freqSamples = freq.sample( 0, length * out.timeToFrame(), out.frameToTime() );
	auto harmonicWeightSamples = harmonicWeights.sample( 0, out.getNumFrames(), out.frameToTime(), 0, out.getNumBins(), out.binToFrequency() );

	std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.getNumFrames() ), [&]( Frame frame )
		{
		const Frequency freq_c = freqSamples[frame];
		if( freq_c <= 1.0 ) return;
		const int numHarmonics = std::floor( out.getHeight() / freq_c );
		for( int harmonic = 0; harmonic < numHarmonics; ++harmonic )
			{
			const float mag = harmonicWeightSamples[buffer_access(harmonic, frame, numHarmonics)] * scale;
			const Frequency harmonicFreq = freq_c * ( harmonic + 1 );
			const Bin bin = harmonicFreq * out.frequencyToBin();
			out.getMF( 0, frame, bin ) = { mag, harmonicFreq };
			}
		} );
	
	return out;
	}

//========================================================================
// Uncategorized
//========================================================================

/*
harmonicFunc must be par_unseq-safe
*/
static PV harmonicScaler( const PV & me, const Func2x1 & series, std::function< Frequency ( Frequency, Harmonic )> harmonicFunc, Harmonic numHarmonics )
	{
	PV out( me.getFormat() );
	out.clearBuffer();

	auto seriesSamples = series.sample( 0, me.getNumFrames(), me.frameToTime(), 0, numHarmonics, 1 );

	for( Channel channel = 0; channel < me.getNumChannels(); ++channel )
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( me.getNumFrames() ), [&]( Frame frame )
			{
			auto seriesSamplesForFrame = seriesSamples.begin() + frame * numHarmonics;

			for( Bin bin = 0; bin < me.getNumBins(); ++bin )
				{
				const MF & source = me.getMF( channel, frame, bin );
				if( source.f <= 1.0f ) continue;
				
				for( Harmonic harmonic = 0; harmonic < numHarmonics; ++harmonic )
					{
					const Frequency harmonicFreq = harmonicFunc( source.f, harmonic );
					const Bin harmonicBin = harmonicFreq * me.frequencyToBin();
					if( harmonicBin >= me.getHeight() ) break;

					MF & dest = out.getMF( channel, frame, harmonicFreq * me.frequencyToBin() );
					const Magnitude overwriteMag = source.m * seriesSamplesForFrame[harmonic];
					if( dest.m < overwriteMag )
						dest = { overwriteMag, harmonicFreq };
					}
				}
			} );

	return out;
	}

PV PV::addOctaves( const Func2x1 & series ) const
	{
	flan_PROCESS_START( PV() );
	return harmonicScaler( *this, series, []( Frequency f, Harmonic h ){ return f * std::pow( 2, h ); }, std::ceil( std::log2( getHeight() ) ) );
	}

PV PV::addHarmonics( const Func2x1 & series ) const
	{
	flan_PROCESS_START( PV() );
	return harmonicScaler( *this, series, []( Frequency f, Harmonic h ){ return f * ( h + 1 ); }, getNumBins() );
	}

PV PV::shape( const Function<MF,MF> & shaper, bool useShiftAlignment ) const
	{
	flan_PROCESS_START( PV() );

	PV out( getFormat() );
	out.clearBuffer();

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		runtimeExecutionPolicyHandler( shaper.getExecutionPolicy(), [&]( auto policy ){
		std::for_each( policy, iota_iter( 0 ), iota_iter( getNumFrames() ), [&]( Frame frame )
			{
			for( Bin bin = 0; bin < getNumBins(); ++bin )
				{
				const MF inMF = getMF( channel, frame, bin );
				const MF shapedMF = shaper( inMF );
				
				if( useShiftAlignment )
					{
					const Bin binShift = bin - inMF.f * frequencyToBin();
					const Bin shapedMFBin = shapedMF.f * frequencyToBin() + binShift;
					if( shapedMFBin < 0 || getNumBins() <= shapedMFBin ) 
						continue;

					MF & outMF = out.getMF( channel, frame, shapedMFBin );
					if( shapedMF.m > outMF.m )
						outMF = shapedMF;
					}
				else 
					{
					out.getMF( channel, frame, bin ) = shapedMF;
					}
				}
			} ); } );
		}
			
	return out;
	}

PV PV::perturb( const Func2x1 & magSigma, const Func2x1 & frqSigma, float damping ) const
	{
	/*
	Velocity base with no const axis is most promising so far.
	Currently being computed with a summation over both x and y of the random sampling for smoothing.
	Other smoothing should give other interesting results.
	Magnitude is currently unimplimented, test other overall process options like simplex and then add mag handling.
	*/

	flan_PROCESS_START( PV() );

	const float epsilon = 0.00001;

	// Function to sample f over domain of this. f must be thread safe or have an appropriate ExecutionPolicy.
	auto sampleOnThis = [&]( const auto & f ){
		return f.sample( 0, getNumFrames(), frameToTime(), 0, getNumBins(), binToFrequency() );
		};

	// Input validation
	auto magSigmaSamples = sampleOnThis( magSigma );
	std::ranges::for_each( magSigmaSamples, []( float & x ){ x = std::max( x, 0.0f ); } );
	auto frqSigmaSamples = sampleOnThis( frqSigma );
	std::ranges::for_each( frqSigmaSamples, []( float & x ){ x = std::max( x, 0.0f ); } );

	// Initialize random engine
	std::default_random_engine rng( std::time( nullptr ) );

	// Sample distributions using random engine. This is done outside the processing loop because random engines aren't thread safe.
	std::vector<float> frqAccels( getNumBins() * getNumFrames() );
	for( Bin i : iota_view( 0u, frqAccels.size() ) )
		{
		const float frqSigma_c = frqSigmaSamples[i];
		frqAccels[i] = frqSigma_c < epsilon ? 0 : std::normal_distribution<float>( 0, frqSigma_c / 20 )( rng );
		}

	std::vector<float> frqVelocs( frqAccels.size() );
	for( Bin bin = 0; bin < getNumBins(); ++bin )
		{
		int prevPos = buffer_access( bin, 0, getNumBins() );
		frqVelocs[prevPos] = frqAccels[prevPos];

		for( Frame frame = 0; frame < getNumFrames(); ++frame )
			{
			const int bufferPos = buffer_access( bin, frame, getNumBins() );
			frqVelocs[bufferPos] = ( frqAccels[bufferPos] + frqVelocs[prevPos] ) * damping;
			prevPos = bufferPos;
			}
		}

	std::vector<float> frqOffsets( frqAccels.size() );
	for( Frame frame = 0; frame < getNumFrames(); ++frame )
		{
		int prevPos = buffer_access( 0, frame, getNumBins() );
		frqOffsets[prevPos] = frqVelocs[prevPos];

		for( Bin bin = 0; bin < getNumBins(); ++bin )
			{
			const int bufferPos = buffer_access( bin, frame, getNumBins() );
			frqOffsets[bufferPos] = ( frqVelocs[bufferPos] + frqOffsets[prevPos] ) * damping;
			prevPos = bufferPos;
			}
		}

	PV out( getFormat() );

	for( Channel channel : iota_view( 0, getNumChannels() ) ) {
		float magOffset  = 0;
		float freqOffset = 0;
		for( Frame frame : iota_view( 0, getNumFrames() ) )	{
			const auto bufferPos = buffer_access( 0, frame, getNumBins() ); 

			const float magSigma_c = magSigmaSamples[bufferPos];
			const float frqSigma_c = frqSigmaSamples[bufferPos];

			magOffset += magSigma_c < epsilon ? 0 : std::normal_distribution<float>( 0, magSigma_c / 20 )( rng );
			//frqOffset += frqSigma_c < epsilon ? 0 : std::normal_distribution<float>( 0, frqSigma_c / 20 )( rng );

			for( Bin bin = 0; bin < getNumBins(); ++bin ) {
				const MF & currentBin = getMF( channel, frame, bin );
				
				
				const Magnitude mag = currentBin.m + magOffset;
				const Frequency freq = currentBin.f + frqOffsets[bin] * 200;
				//const Frequency freq = currentBin.f * exp2( frqOffsets[frame] );
				out.getMF( channel, frame, bin ) = { mag, freq };
				}
			}
		}

	return out;
	}

PV predicateNLoudestPartials( const PV & me, const Function<Time, Bin> & numBins, std::function< bool (Bin, Bin) > predicate )
	{
	flan_FUNCTION_LOG;

	// Input validation
	auto numBinsSamples = numBins.sample( 0, me.getNumFrames(), me.frameToTime() );
	std::ranges::for_each( numBinsSamples, [&]( Bin & b ){ b = std::clamp( b, 0, me.getNumFrames() ); } );

	PV out( me.getFormat() );
	out.clearBuffer();

	for( Channel channel = 0; channel < me.getNumChannels(); ++channel )
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( me.getNumFrames() ), [&]( Frame frame )
			{
			// Create a buffer storing the frame's magnitudes and which bin that magnitude came from
			using BinMag = std::pair<Bin, Magnitude>;
			std::vector<BinMag> indexAndVolumes( me.getNumBins() );
			for( Bin bin = 0; bin < me.getNumBins(); ++bin )
				indexAndVolumes[bin] = { bin, me.getMF( channel, frame, bin ).m };

			// Sort that buffer by magnitudes
			std::sort( std::execution::par_unseq, indexAndVolumes.begin(), indexAndVolumes.end(), []( BinMag & a, BinMag & b )
				{ return abs(a.second) > abs(b.second); } );
				
			// The sorted buffer now contains bin indices in loudness order, so go back through and use the predicate 
			// to decide which to copy to the output.
			for( Bin bin = 0; bin < me.getNumBins(); ++bin )
				{
				const Bin actualBin = indexAndVolumes[bin].first;
				const MF currentMF = me.getMF( channel, frame, actualBin );
				if( predicate( bin, numBinsSamples[frame] ) )
					out.setMF( channel, frame, actualBin, currentMF );
				else
					out.setMF( channel, frame, actualBin, { 0.0, currentMF.f } );
				}
			} );

	return out;
	}

PV PV::retainNLoudestPartials( const Function<Time, Bin> & numBins ) const
	{
	flan_PROCESS_START( PV() );
	return predicateNLoudestPartials( *this, numBins, []( Bin a, Bin b ){ return a < b; } );
	}

PV PV::removeNLoudestPartials( const Function<Time, Bin> & numBins ) const
	{
	flan_PROCESS_START( PV() );
	return predicateNLoudestPartials( *this, numBins, []( Bin a, Bin b ){ return a >= b; } );
	}

PV PV::resonate( Time length, const Func2x1 & decay ) const
	{
	flan_PROCESS_START( PV() );

	// Input validation
	if( length < 0 ) 
		length = 0;

	auto decaySamples = decay.sample( 0, getNumFrames(), frameToTime(), 0, getNumBins(), binToFrequency() );
	std::ranges::for_each( decaySamples, []( float & x ){ x = std::max( x, 0.0f ); } );

	auto format = getFormat();
	format.numFrames = getNumFrames() + ceil( timeToFrame() * length );
	PV out( format );

	// Copy first frame into out
	for( Channel channel = 0; channel < out.getNumChannels(); ++channel )
		for( Bin bin = 0; bin < out.getNumBins(); ++bin )
			{
			out.setMF( channel, 0, bin, getMF( channel, 0, bin ) );
			}
	
	const float secondsPerFrame_c = frameToTime();

	for( Channel channel = 0; channel < out.getNumChannels(); ++channel )
		std::for_each( std::execution::par_unseq, iota_iter( 1 ), iota_iter( out.getNumFrames() ), [&]( Frame frame )
			{
			for( Bin bin = 0; bin < out.getNumBins(); ++bin )
				{
				const float currentDecay = std::pow( decaySamples[buffer_access( bin, frame, getNumBins() )], secondsPerFrame_c );
				const float decayedAmp = out.getMF( channel, frame - 1, bin ).m * currentDecay;
				if( frame < getNumFrames() && getMF( channel, frame, bin ).m > decayedAmp )
					out.setMF( channel, frame, bin, getMF( channel, frame, bin ) );
				else
					out.setMF( channel, frame, bin, { decayedAmp, out.getMF( channel, frame - 1, bin ).f } );
				}
			} );

	return out;
	}

} //End namespace flan
