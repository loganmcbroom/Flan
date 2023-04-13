#include "flan/Audio/Audio.h"

#include <iostream>
#include <algorithm>
#include <complex>
#include <numeric>
#include <execution> 
#include <ranges>

#include "flan/WindowFunctions.h"
#include "WDL/resample.h"
#include "flan/Utility/iota_iter.h"
#include "flan/Utility/execution.h"

using namespace flan;
using namespace std::ranges;

Audio Audio::copy() const { return AudioBuffer::copy(); }

Audio::Audio() 
	: AudioBuffer() 
	{}

Audio::Audio( std::vector<float> && buffer, Channel numChannels, SampleRate sr )
	: AudioBuffer( std::move( buffer ), numChannels, sr )
	{}

Audio::Audio( const AudioBuffer::Format & other ) 
	: AudioBuffer( other ) 
	{}

Audio::Audio( const std::string & filename ) 
	: AudioBuffer( filename ) 
	{}

Audio::Audio( AudioBuffer && other )
	: AudioBuffer( std::move( other ) )
	{
	}

//========================================================
// Helper Functions
//========================================================

// static bool doSampleRatesMatch( const std::vector<Audio> & ins )
// 	{
// 	//Check if all sample rates match the first file
// 	if( ins.empty() ) return true;
// 	const uint32_t sampleRate = ins[0].getSampleRate();
// 	for( const Audio & in : ins ) 
// 		if( in.getSampleRate() != sampleRate )
// 			{
// 			std::cout << "Mismatched sample rates" << std::endl;
// 			return false;
// 			}
// 	return true;
// 	}

static bool doChannelCountsMatch( const std::vector<Audio> & ins )
	{
	//Check if all channel counts match the first one
	if( ins.empty() ) return true;
	const Channel numchannels = ins[0].getNumChannels();
	for( auto & in : ins ) 
		if( in.getNumChannels() != numchannels )
			return false;
	return true;
	}

static Channel getMaxNumChannels( const std::vector<const Audio *> & ins )
	{
	const auto maxNumChannelsIter = std::max_element( ins.begin(), ins.end(), []( const Audio * a, const Audio * b )
		{ 
		return a->getNumChannels() < b->getNumChannels();
		} );

	return (*maxNumChannelsIter)->getNumChannels();
	}

static Frame getMaxNumFrames( const std::vector<Audio> & ins )
	{
	return std::max_element( ins.begin(), ins.end(), []( const Audio & a, const Audio & b )
		{ 
		return a.getNumFrames() < b.getNumFrames();
		} )->getNumFrames();
	}

template<typename T>
static std::vector<const T *> getPtrs( const std::vector<T> & ts )
	{
	std::vector<const T *> ptrs( ts.size() );
	std::transform( ts.begin(), ts.end(), ptrs.begin(), []( const T & t ){ return &t; } );
	return std::move( ptrs );
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

Audio Audio::invertPhase() const
	{
	flan_PROCESS_START( Audio() );
	Audio out = copy();
	std::for_each( std::execution::par_unseq, out.getBuffer().begin(), out.getBuffer().end(), []( Sample & s ){ s = -s; } );
	return out;
	}

Audio Audio::modifyVolume( const Func1x1 & volumeLevel ) const
	{
	flan_PROCESS_START( Audio() );
	Audio out = copy();
	out.modifyVolumeInPlace( volumeLevel );
	return out;
	}

void Audio::modifyVolumeInPlace( const Func1x1 & volumeLevel )
	{
	flan_PROCESS_START();

	const auto sampledVolumeLevel = volumeLevel.sample( 0, getNumFrames(), frameToTime( 1 ) );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( getNumFrames() ), [&]( Frame frame )
			{
			getSample( channel, frame ) *= sampledVolumeLevel[frame];
			} );
	}

Audio Audio::setVolume( const Func1x1 & level ) const
	{
	flan_PROCESS_START( Audio() );

	// Divide by getMaxSampleMagnitude to normalize, multiply by level to set
	const Sample maxMag = getMaxSampleMagnitude();
	if( maxMag == 0 ) return copy();
	return modifyVolume( [&]( float t ){ return level(t) / maxMag; } );
	}

Audio Audio::shift( Time shift ) const
	{
	flan_PROCESS_START( Audio() );

	if( shift <= 0 ) return cut( -shift, getLength() );

	const Frame shiftFrames = timeToFrame( shift );

	auto format = getFormat();
	format.numFrames += shiftFrames;
	Audio out( format );
	out.clearBuffer();

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( getNumFrames() ), [&]( Frame frame )
			{
			out.setSample( channel, frame + shiftFrames, getSample( channel, frame ) );
			} );

	return out;
	}

Audio Audio::waveshape( const Func1x1 & shaper ) const
	{
	flan_PROCESS_START( Audio() );

	Audio out( getFormat() );
	std::transform( std::execution::par_unseq, getBuffer().begin(), getBuffer().end(), out.getBuffer().begin(), 
		[&]( Sample s ){ return shaper(s); } ); // Can't just pass shaper since transform will take it by value
	return out;
	}

// Audio Audio::pan( Function<Time,vec2> panPosition, Function<Time, vec2> listenerPosition, std::vector<vec2> speakerPositions, bool haasEffect ) const
// 	{
	/*
	Panning is a much more complex process than anyone gives it credit for. Panning two speakers isn't difficult, but
	once you have three there are an infinite number of ways to weigh them to produce a sound appearing to come from a 
	particular position. In linear algebra terms, you have a basis of three (or more) vectors for a 2d space. Selecting a 
	single weighting requires an additional constraint. 

	This currently chooses no additional constraint as it is unimplimented, but a number of constraints are discussed here:
	https://doc.flux.audio/en_US/spat_revolution_doc/Spatialisation_Technology_Panning_Algorithms.html
	*/

// 	flan_PROCESS_START( Audio() );

// 	Audio out = copy();

// 	const std::vector<vec2> panPositionSamples = panPosition.sample( 0, getNumFrames(), frameToTime() );
// 	const std::vector<vec2> listenerPositionSamples = listenerPosition.sample( 0, getNumFrames(), frameToTime() );
	
// 	for( Channel channel = 0; channel < getNumChannels(); ++channel )
// 		{
// 		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( getNumFrames() ), [&]( Frame frame )
// 			{
// 			// Subtracting the current listener position from all the other vec2s is equivalent to moving the listener
// 			const vec2 relativePanPos = panPositionSamples[frame] - listenerPositionSamples[frame];
// 			const vec2 relativeSpeakerPos = speakerPositions[channel] - listenerPositionSamples[frame];

			
// 			} );
// 		}

// 	return out;
// 	}

Audio Audio::pan( const Func1x1 & panAmount ) const
	{
	flan_PROCESS_START( Audio() );

	if( getNumChannels() != 1 && getNumChannels() != 2 )
		return Audio();

	Audio out = getNumChannels() == 1 ? convertToStereo() : copy();
	out.panInPlace( panAmount );
	return out;
	}

void Audio::panInPlace( const Func1x1 & panAmount )
	{
	flan_PROCESS_START();

	if( getNumChannels() != 2 ) return;

	const auto sampledPanAmount = panAmount.sample( 0, getNumFrames(), frameToTime( 1 ) );

	for( Channel channel = 0; channel < 2; ++channel )
		{
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( getNumFrames() ), [&]( Frame frame )
			{
			const float pan = sampledPanAmount[frame] / 2.0f + 1.0f; // Convert [-1,1] to [0,1]
			const float panFuncInput = channel == 0 ? pan : ( 1.0f - pan );
			const float scale = Interpolators::sine2( panFuncInput );
			getSample( channel, frame ) *= scale;
			} );
		}
	}

Audio Audio::widen( const Func1x1 & widenAmount ) const
	{
	flan_FUNCTION_LOG
	return convertToMidSide().pan( widenAmount ).convertToLeftRight();
	}

Audio Audio::reverse() const
	{
	flan_PROCESS_START( Audio() );

	Audio out( getFormat() );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		const auto iter = getBuffer().begin() + channel * getNumFrames();
		std::copy( std::execution::par_unseq, iter, iter + getNumFrames(), out.getBuffer().rbegin() + channel * getNumFrames() );
		}

	return out;
	}

Audio Audio::cut( Time startTime, Time endTime, Time startFade, Time endFade ) const
	{
	flan_PROCESS_START( Audio() );

	return cutFrames( 
		timeToFrame( startTime 	), 
		timeToFrame( endTime 	), 
		timeToFrame( startFade 	), 
		timeToFrame( endFade 	)
		);
	}

Audio Audio::cutFrames( Frame start, Frame end, Frame startFade, Frame endFade ) const
	{
	flan_PROCESS_START( Audio() );

	// Input validation
	if( end <= start ) return Audio();
	const Frame startFrame = std::clamp( start, 0, getNumFrames() - 1 );
	const Frame endFrame   = std::clamp( end,   0, getNumFrames() - 1 );
	// Any fade errors are validated in Audio::fades
	
	auto format = getFormat();
	format.numFrames = end - start;
	Audio out( format );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.getNumFrames() ), [&]( Frame frame )
			{
			out.setSample( channel, frame, getSample( channel, start + frame ) );
			} );

	out.fadeFramesInPlace( startFade, endFade, Interpolators::sqrt );
	return out;
	}

Audio Audio::repitch( const Func1x1 & factor, Time granulTime, WDLResampleType quality ) const
	{
	flan_PROCESS_START( Audio() );

	// Input validation
	Frame granul = timeToFrame( granulTime );
	if( granul < 1 ) granul = 1;

	// This will be inverted momentarily, hence the name
	auto factorSamplesInv = factor.sample( 0, std::ceil( getNumFrames() / float( granul ) ), granulTime );
	std::for_each( std::execution::par_unseq, factorSamplesInv.begin(), factorSamplesInv.end(), []( float & sample ){
		constexpr float bound = 1000.0f; // Arbitrary
		sample = std::clamp( 1.0f / sample, 1.0f / bound, bound ); 
		} );

	// Estimate output length
	const Frame numOutFrames = std::ceil( std::accumulate( factorSamplesInv.begin(), factorSamplesInv.end(), 0.0f ) * granul );

	auto format = getFormat();
	format.numFrames = numOutFrames;
	Audio out( format );
	
	WDL_Resampler rs;
		 if( quality == WDLResampleType::Sinc ) 			rs.SetMode( true,  0, true, 64 ); // reasonable sinc resampling
	else if( quality == WDLResampleType::Linear ) 			rs.SetMode( true,  1, false    ); // linear interpolation with simple filter
	else if( quality == WDLResampleType::Uninterpolated ) 	rs.SetMode( false, 0, false    ); // no interpolation, useful for dirty sound

	const Channel numChannels = getNumChannels();
	const Frame inFrames = getNumFrames();
	std::vector<Sample> rsoutbuf( numChannels * granul );

	Frame inFrame = 0;
	Frame outFrame = 0;
	while( inFrame < getNumFrames() )
		{
		const double factor_to_use = factorSamplesInv[ std::floor( inFrame / float( granul ) ) ];
		rs.SetRates( getSampleRate(), double( getSampleRate() ) * factor_to_use );
		WDL_ResampleSample* rsinbuf = nullptr;
		const Frame wanted = rs.ResamplePrepare( granul, numChannels, &rsinbuf );

		for( Channel channel = 0; channel < numChannels; ++channel )
			for( Frame frame = 0; frame < wanted; ++frame )
				{
				if( inFrame + frame < inFrames )
					rsinbuf[ frame * numChannels + channel ] = getSample( channel, inFrame + frame );
				else
					rsinbuf[ frame * numChannels + channel ] = 0.0;
				}

		// Process
		rs.ResampleOut( rsoutbuf.data(), wanted, granul, numChannels );

		// Copy resampled data to output buffer
		for( Channel channel = 0; channel < numChannels; ++channel )
			for( Frame frame = 0; frame < granul; ++frame )
				if( outFrame + frame < numOutFrames )
					out.setSample( channel, outFrame + frame, rsoutbuf[ frame * numChannels + channel ] );

		outFrame += granul;
		inFrame += wanted;
		}

	return out;
	}

Audio Audio::convolve( const Function<Time, std::vector<float>> & ir ) const
	{
	flan_PROCESS_START( Audio() );

	auto irSamples = ir.sample( 0, getNumFrames(), frameToTime( 1 ) );
	const size_t maxIrSize = max_element( irSamples, less(), []( const std::vector<float> & v ){ return v.size(); } )->size();

	auto format = getFormat();
	format.numFrames = getNumFrames() + maxIrSize; // This may be more frames than needed if the ir size changes
	Audio out( format );

	for( Channel channel = 0; channel < out.getNumChannels(); ++channel )
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.getNumFrames() ), [&]( Frame outFrame )
			{
			const auto ir_c = irSamples[outFrame];

			Sample convolutionSum = 0;
			for( Frame irFrame = 0; irFrame < ir_c.size(); ++irFrame )
				{
				const Frame inFrame = outFrame + 1 - irFrame;
				if( 0 <= inFrame && inFrame < getNumFrames() )
					convolutionSum += ir_c[irFrame] * getSample( channel, inFrame );
				}

			out.setSample( channel, outFrame, convolutionSum );
			} );

	return out;
	}

Audio Audio::convolve( const Audio & ir ) const
	{
	flan_PROCESS_START( Audio() );
	if( ir.isNull() ) return Audio();

	auto format = getFormat();
	format.numFrames = getNumFrames() + ir.getNumFrames();
	Audio out( format );

	for( Channel channel = 0; channel < out.getNumChannels(); ++channel )
		{
		const Channel irChannel = channel % ir.getNumChannels();
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.getNumFrames() ), [&]( Frame outFrame )
			{
			Sample convolutionSum = 0;
			for( Frame irFrame = 0; irFrame < ir.getNumFrames(); ++irFrame )
				{
				const Frame inFrame = outFrame + 1 - irFrame;
				if( 0 <= inFrame && inFrame < getNumFrames() )
					convolutionSum += ir.getSample( irChannel, irFrame ) * getSample( channel, inFrame );
				}

			out.setSample( channel, outFrame, convolutionSum );
			} );
		}

	return out;
	}

Audio Audio::iterate( uint32_t n, const AudioMod & mod, bool feedback ) const
	{
	flan_PROCESS_START( Audio() );
	// Half a length is removed to avoid length rounding changing the number of events generated
	return texture( getLength() * ( float( n ) - 0.5f ), 1.0f / getLength(), 0, mod, feedback );
	}

Audio Audio::delay( Time length, const Func1x1 & delayTime, const Func1x1 & decay, const AudioMod & mod ) const
	{
	flan_PROCESS_START( Audio() );

	if( length <= 0 ) return Audio();

	return Audio();

	auto delayTimeSamples = delayTime.sample( 0, timeToFrame( length ), frameToTime( 1 ) );
	auto decaySamples = decay.sample( 0, timeToFrame( length ), frameToTime( 1 ) );

	auto eventsPerSecond = [&]( float t )
		{ 
		const float delayTime_c = delayTimeSamples[ timeToFrame( t ) ];
		if( delayTime_c <= 0 ) return 1.0f / float( getSampleRate() );
		return 1.0f / delayTime_c;
		};

	const AudioMod delayMod = [&]( const Audio & in, float t ) -> Audio
		{ 
		if( t == 0 ) return in.copy();

		Audio out = mod ? mod( in, t ) : in.copy();

		// Modify gain
		const float currentDecay = decaySamples[ timeToFrame( t ) ];
		std::for_each( out.getBuffer().begin(), out.getBuffer().end(), [currentDecay]( Sample & s ){ s *= currentDecay; } );

		return out;
		};

	return texture( length, eventsPerSecond, 0, delayMod, true );
	}

Audio Audio::fade( Time start, Time end, Interpolator interp ) const
	{
	flan_PROCESS_START( Audio() );
	return fadeFrames( timeToFrame( start ), timeToFrame( end ), interp );
	}

void Audio::fadeInPlace( Time start, Time end, Interpolator interp )
	{
	flan_PROCESS_START();
	fadeFramesInPlace( timeToFrame( start ), timeToFrame( end ), interp );
	}

Audio Audio::fadeFrames( Frame start, Frame end, Interpolator interp ) const
	{
	flan_PROCESS_START( Audio() );
	Audio out = copy();
	out.fadeFramesInPlace( start, end, interp );
	return out;
	}

void Audio::fadeFramesInPlace( Frame start, Frame end, Interpolator interp )
	{
	flan_PROCESS_START();

	// Input validation
	start = std::max( 0, start );
	end   = std::max( 0, end   );
	if( start + end > getNumFrames() )
		{
		const float scale = float( getNumFrames() ) / ( start + end );
		start = std::floor( start * scale );
		end	= std::floor( end * scale );
		}
	if( start == 0 && end == 0 ) return;

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		for( Frame frame = 0; frame < start; ++frame )
			{
			const float fadeAmt = interp( float( frame ) / start );
			getSample( channel, frame ) *= fadeAmt;
			}
		for( Frame frame = 0; frame < end; ++frame )
			{
			const float fadeAmt = interp( float( frame ) / end );
			getSample( channel, getNumFrames() - 1 - frame ) *= fadeAmt;
			}
		}
	}

// Audio Audio::lowPass( Func1x1 cutoff, uint32_t taps, flan_CANCEL_ARG_CPP ) const
// 	{
// 	flan_PROCESS_START( Audio() );

// 	if( taps == 0 ) return copy();

// 	std::vector<Func1x1> ir;
// 	for( uint32_t tap = 0; tap <= taps; ++tap ) 
// 		ir.push_back( [&cutoff, this, taps, tap]( float t )
// 			{
// 			const float n = float(tap) - int(taps) / 2.0f;
// 			const float N = getSampleRate();
// 			const float K = 2.0f * cutoff(t);
// 			if( n == 0 ) return K / N;
// 			else return (float) ( std::sin( pi * n * K / N ) / ( N * sin( pi * n / N ) ) * Windows::Hann( float(tap) / float(taps) ) ); 
// 			} );

// 	return convolve( ir, canceller );
// 	}

Audio Audio::grainSelect( Time length, const Func1x1 & grainsPerSecond, const Func1x1 & scatter, 
		const Func1x1 & selection, const Func1x1 & grainLength, const Func1x1 & fade, const AudioMod & mod ) const
	{
	flan_PROCESS_START( Audio() );

	// Input validation
	if( length <= 0 ) return Audio();

	const auto selectionSamples = selection		.sample( 0, timeToFrame( length ), frameToTime( 1 ) );
	const auto grainLengthSamples = grainLength	.sample( 0, timeToFrame( length ), frameToTime( 1 ) );
	const auto fadeSamples 		= fade			.sample( 0, timeToFrame( length ), frameToTime( 1 ) );

	AudioMod grainMod = [&]( const Audio & in, float t ) -> Audio
		{
		const Time selection_c = selectionSamples[ timeToFrame( t ) ];
		const Time grainLength_c = grainLengthSamples[ timeToFrame( t ) ];
		const Time fade_c = fadeSamples[ timeToFrame( t ) ];
		Audio grain = in.cut( selection_c, selection_c + grainLength_c, fade_c, fade_c );
		return std::move( grain );
		};

	return texture( length, grainsPerSecond, scatter, mod ? mod( grainMod ) : grainMod, false );
	}

std::vector<Audio> Audio::chop( Time sliceLength, Time fade ) const
	{
	flan_PROCESS_START( std::vector<Audio>() );

	// Input validation
	if( sliceLength <= 0 ) return std::vector<Audio>();

	// Get all but last slice
	std::vector<Audio> outs;
	Time slicePosition = 0;
	while( slicePosition + sliceLength < getLength() )
		{
		outs.push_back( cut( slicePosition, slicePosition + sliceLength, fade, fade ) );
		slicePosition += sliceLength;
		}

	return outs;
	}

Audio Audio::rearrange( Time sliceLength, Time fade ) const
	{
	flan_PROCESS_START( Audio() );

	// + fade here accounts for + fade / 2 at both ends
	std::vector<Audio> chops = chop( sliceLength + fade, fade );
	if( chops.size() < 2 )
		return Audio();

	chops.pop_back(); // Final slice usually isn't the correct length

	std::random_device rd;
	std::mt19937 g( rd() );
	std::shuffle( chops.begin(), chops.end(), g );

	return Audio::join( std::move( chops ), -fade );
	}
	

//========================================================
// Multi-In Procs
//========================================================

Audio Audio::mix( std::vector<const Audio *> ins, std::vector<Time> startTimes, const std::vector<const Func1x1 *> & gains )
	{
	flan_FUNCTION_LOG;

	// Input validation
	if( ins.empty() ) return Audio();

	// If startTimes is too small fill with 0
	for( auto _ : iota_view( startTimes.size(), ins.size() ) )
		startTimes.push_back( 0 );

	// Convert start times to frames
	std::vector<Frame> startFrames;
	transform( startTimes, std::back_inserter( startFrames ), [&]( Time t ) -> Frame { 
		return ins[0]->timeToFrame( t ); 
		} );

	// Sample each function, if there aren't enough use constant 1
	// We only sample the gains when they are needed, but because these functions are relative to global time 0, we need to be 
	// careful later when accessing them
	std::vector<std::vector<float>> gainsSamples;
	const Func1x1 oneFunc = 1;
	std::transform( iota_iter( 0 ), iota_iter( ins.size() ), std::back_inserter( gainsSamples ), [&]( Index i ) {
		const Func1x1 * f = i < gains.size() ? gains[i] : &oneFunc;
		return f->sample( startFrames[i], startFrames[i] + ins[i]->getNumFrames(), ins[i]->frameToTime( 1 ) );
		} );

	// Output setup
	auto format = ins[0]->getFormat();
	// Get the maximum number of channels among inputs
	format.numChannels = (*max_element( ins, less(), []( const Audio * p ) { return p->getNumChannels(); } ) )->getNumChannels();
	// Get the maximum of output sizes required by each input
	format.numFrames = 0;
	for( const Index i : iota_view( 0u, ins.size() ) ) {
		const Frame framesNeededForI = ins[i]->getNumFrames() + startFrames[i];
		format.numFrames = std::max( format.numFrames, framesNeededForI );
		}
	Audio out( format );
	out.clearBuffer();

	// If any inputs don't match the first sample rate, they need to be converted
	// This vector just gives the resampled inputs somewhere to live until they get mixed and can be deleted
	std::vector<Audio> resampledInputs;
	// I'm not going to claim this pointer reference is indicative of good design, but it's a consequence of c++ not handling reference containers
	for( const Audio *& a : ins ) 
		if( a->getSampleRate() != out.getSampleRate() )
			{
			resampledInputs.push_back( a->resample( out.getSampleRate() ) );
			a = &resampledInputs.back();
			}

	// For each input, write to output
	for( Index in = 0; in < ins.size(); ++in )
		{
		const Audio & me = *ins[in];
		for( Channel channel = 0; channel < me.getNumChannels(); ++channel )
			std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( me.getNumFrames() ), [&]( Frame localFrame )
				{
				const Frame globalFrame = localFrame + startFrames[in];
				if( 0 <= globalFrame && globalFrame < format.numFrames ) {
					const float gain = gainsSamples[in][localFrame];
					out.getSample( channel, globalFrame ) += me.getSample( channel, localFrame ) * gain;
					}
				} );
		}

	return out;
	}

Audio Audio::mix( const std::vector<Audio> & ins, const std::vector<Time> & startTimes, const std::vector<Func1x1> & gains )
	{
	return mix( getPtrs( ins ), startTimes, getPtrs( gains ) );
	}

Audio Audio::join( const std::vector<const Audio *> & ins, Time offset )
	{
	flan_FUNCTION_LOG;

	if( ins.empty() ) return Audio();

	// Get input Audio lengths
	std::vector<Time> jumps( { 0 } );
	transform( ins, std::back_inserter( jumps ), []( const Audio * a ){ return a->getLength(); } );

	// Sum jumps to get mix positions
	std::vector<Time> startTimes;
	std::partial_sum( jumps.begin(), jumps.end() - 1, std::back_inserter( startTimes ), [&]( Time a, Time b ){ 
		return a + b + offset; 
		} );

	return mix( ins, startTimes, {} );
	}

Audio Audio::join( const std::vector<Audio> & ins, Time offset )
	{
	return join( getPtrs( ins ), offset );
	}

Audio Audio::select( const std::vector<const Audio *> & ins, const Func1x1 & selection, const std::vector<Time> & startTimes )
	{
	flan_FUNCTION_LOG;

	// Generate balances from selection
	std::vector<Func1x1> balances;
	for( Index i = 0; i < ins.size(); ++i )
		{
		balances.emplace_back( [&selection, i]( Time t )
			{
			const float distance = std::abs( selection( t ) - i );
			if( distance >= 1 ) return 0.0f;
			else return std::sqrt( 1.0f - distance );
			} );
		}
		
	return mix( ins, startTimes, getPtrs( balances ) );
	}

Audio Audio::select( const std::vector<Audio> & ins, const Func1x1 & selection, std::vector<Time> startTimes )
	{
	return select( getPtrs( ins ), selection, startTimes );
	}