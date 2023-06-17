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
#include "Audio.h"

using namespace flan;
using namespace std::ranges;

Audio Audio::copy() const { return AudioBuffer::copy(); }

Audio::Audio() 
	: AudioBuffer() 
	{}

Audio::Audio( std::vector<float> && buffer, Channel numChannels, FrameRate sr )
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

// std::vector<Audio> Audio::splitChannels() const
// 	{
// 	auto format = getFormat();
// 	format.numChannels = 1;

// 	std::vector<Audio> out( getNumChannels(), format );

// 	for( Channel channel = 0; channel < getNumChannels(); ++channel )
// 		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( getNumFrames() ), [&]( Frame frame )
// 			{
// 			out[channel].getSample( 0, frame ) = getSample( channel, frame );
// 			} );

// 	return out;
// 	}

Audio Audio::combineChannels( const Audio & other ) const
	{
	Audio::Format format;
	format.numChannels = getNumChannels() + other.getNumChannels();
	format.numFrames = std::max( getNumFrames(), other.getNumFrames() );
	format.sampleRate = getSampleRate();
	Audio out( format );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( getNumFrames() ), [&]( Frame frame )
			{ out.getSample( channel, frame ) = getSample( channel, frame ); } );

	for( Channel channel = 0; channel < other.getNumChannels(); ++channel )
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( other.getNumFrames() ), [&]( Frame frame )
			{ out.getSample( channel + getNumChannels(), frame ) = other.getSample( channel, frame ); } );
		
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

Audio Audio::shift( Second shift ) const
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

// Audio Audio::pan( Function<Second,vec2> panPosition, Function<Second, vec2> listenerPosition, std::vector<vec2> speakerPositions, bool haasEffect ) const
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

Audio Audio::cut( Second startTime, Second endTime, Second startFade, Second endFade ) const
	{
	flan_PROCESS_START( Audio() );

	return cutFrames( 
		timeToFrame( startTime	), 
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

Audio Audio::repitch( const Func1x1 & factor, Second granulTime, WDLResampleType quality ) const
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

Audio Audio::convolve( const Function<Second, std::vector<float>> & ir ) const
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

Audio Audio::delay( Second length, const Func1x1 & delayTime, const Func1x1 & decay, const AudioMod & mod ) const
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

Audio Audio::fade( Second start, Second end, Interpolator interp ) const
	{
	flan_PROCESS_START( Audio() );
	return fadeFrames( timeToFrame( start ), timeToFrame( end ), interp );
	}

void Audio::fadeInPlace( Second start, Second end, Interpolator interp )
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

Audio Audio::grainSelect( Second length, const Func1x1 & grainsPerSecond, const Func1x1 & scatter, 
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
		const Second selection_c = selectionSamples[ timeToFrame( t ) ];
		const Second grainLength_c = grainLengthSamples[ timeToFrame( t ) ];
		const Second fade_c = fadeSamples[ timeToFrame( t ) ];
		Audio grain = in.cut( selection_c, selection_c + grainLength_c, fade_c, fade_c );
		return std::move( grain );
		};

	return texture( length, grainsPerSecond, scatter, mod ? mod( grainMod ) : grainMod, false );
	}

std::vector<Audio> Audio::chop( Second sliceLength, Second fade ) const
	{
	flan_PROCESS_START( std::vector<Audio>() );

	// Input validation
	if( sliceLength <= 0 ) return std::vector<Audio>();

	// Get all but last slice
	std::vector<Audio> outs;
	Second slicePosition = 0;
	while( slicePosition + sliceLength < getLength() )
		{
		outs.push_back( cut( slicePosition, slicePosition + sliceLength, fade, fade ) );
		slicePosition += sliceLength;
		}

	return outs;
	}

Audio Audio::rearrange( Second sliceLength, Second fade ) const
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

Audio Audio::compress( 
	Function<Second, Decibel> threshold, 
	Function<Second, float> ratio, 
	Function<Second, Second> attack, 
	Function<Second, Second> release, 
	Function<Second, Decibel> knee_width, 
	const Audio * sidechain_source ) const
	{
	/*
	For details on compression, both in general and as referenced in this function:
	"Digital Dynamic Range Compressor Design — A Tutorial and Analysis"
	https://www.eecs.qmul.ac.uk/~josh/documents/2012/GiannoulisMassbergReiss-dynamicrangecompression-JAES2012.pdf
	*/

	flan_PROCESS_START( Audio() );

	if( sidechain_source == nullptr )
		sidechain_source = this;

	Audio out = copy();

	// Sample all inputs
	auto threshold_sampled 	= threshold .sample( 0, getNumFrames(), frameToTime( 1 ) );
	auto ratio_sampled 		= ratio     .sample( 0, getNumFrames(), frameToTime( 1 ) );
	auto attack_sampled 	= attack    .sample( 0, getNumFrames(), frameToTime( 1 ) );
	auto release_sampled 	= release   .sample( 0, getNumFrames(), frameToTime( 1 ) );
	auto knee_width_sampled = knee_width.sample( 0, getNumFrames(), frameToTime( 1 ) );

	// (4)
	auto gain_computer = [&]( Decibel x_G, Decibel threshold, Decibel knee_width, float ratio ) -> Decibel
		{
		const Decibel overshoot = x_G - threshold;
		if( overshoot <= -knee_width/2.0f )
			return x_G;
		else if( overshoot >= knee_width / 2.0f )
			return x_G + overshoot * ( 1 / ratio - 1 );
		else
			{
			const Decibel z = overshoot + knee_width/2.0f;
			return x_G + ( 1 / ratio - 1 ) * z*z / ( 2.0f * knee_width );
			}
		};

	// (7)
	auto time_to_alpha = [sr = getSampleRate()]( Second t ){ return std::exp( -1.0f / ( t * sr ) ); };

	// (17)
	// The term "peak detector" is used in the literature. It is more or less a filter.
	auto smooth_decouple_peak_detector = [&]( float & y_1, float & y_L, Frame f, float x_L )
		{
		float a_A = time_to_alpha( attack_sampled[ f ] );
		float a_R = time_to_alpha( release_sampled[ f ] );
		y_1 = std::max( x_L, a_R * y_1 );
		y_L = a_A * y_L + ( 1.0f - a_A ) * y_1;
		return y_L;
		};

	std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( getNumChannels() ), [&]( Channel channel)
		{
		// Peak detector buffers
		float y_1 = 0;
		float y_L = 0;
		for( Frame frame = 0; frame < getNumFrames(); ++frame )
			{
			// (23)
			const Sample x = sidechain_source->getSample( channel, frame );
			const Decibel x_G = 20.0f * std::log10( std::max( std::abs( x ), 1e-6f ) ); // Max is to avoid -inf from log
			const Decibel y_G = gain_computer( x_G, 
				threshold_sampled[ frame ], 
				knee_width_sampled[ frame ], 
				ratio_sampled[ frame ] );
			const Decibel x_L = x_G - y_G;
			const Decibel c_dB = -smooth_decouple_peak_detector( y_1, y_L, frame, x_L );

			const Sample c = std::pow( 10.0f, c_dB / 20.0f );
			out.getSample( channel, frame ) *= c;
			}
		} );

	return out;
	}

// Audio Audio::limit(
// 	Function<Second, Decibel> threshold, 
// 	Function<Second, float> compression_ratio, 
// 	Function<Second, Second> attack, 
// 	Function<Second, Second> release, 
// 	Function<Second, Decibel> knee_width ) const
// 	{
// 	return Audio();
// 	}

//========================================================
// Multi-In Procs
//========================================================

Audio Audio::mix( 
	std::vector<const Audio *> ins, 
	std::vector<Second> startTimes, 
	const std::vector<const Function<Second, Amplitude> *> & amplitudes 
	)
	{
	flan_FUNCTION_LOG;

	// Input validation
	if( ins.empty() ) return Audio();

	// If startTimes is too small fill with 0
	for( auto _ : iota_view( startTimes.size(), ins.size() ) )
		startTimes.push_back( 0 );

	// Convert start times to frames
	std::vector<Frame> startFrames;
	transform( startTimes, std::back_inserter( startFrames ), [&]( Second t ) -> Frame { 
		return ins[0]->timeToFrame( t ); 
		} );

	// Sample each function, if there aren't enough use constant 1
	// We only sample the gains when they are needed, but because these functions are relative to global time 0, we need to be 
	// careful later when accessing them
	std::vector<std::vector<Amplitude>> amplitude_samples;
	const Function<Second, Amplitude> level_func = 1;
	std::transform( iota_iter( 0 ), iota_iter( ins.size() ), std::back_inserter( amplitude_samples ), [&]( Index i ) {
		const Function<Second, Amplitude> * f = i < amplitudes.size() ? amplitudes[i] : &level_func;
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
					const float gain = amplitude_samples[in][localFrame];
					out.getSample( channel, globalFrame ) += me.getSample( channel, localFrame ) * gain;
					}
				} );
		}

	return out;
	}

Audio Audio::mix( 
	const std::vector<Audio> & ins, 
	const std::vector<Second> & startTimes, 
	const std::vector<Function<Second, Amplitude>> & amplitudes 
	)
	{
	return mix( getPtrs( ins ), startTimes, getPtrs( amplitudes ) );
	}

Audio Audio::join( const std::vector<const Audio *> & ins, Second offset )
	{
	flan_FUNCTION_LOG;

	if( ins.empty() ) return Audio();

	// Get input Audio lengths
	std::vector<Second> jumps( { 0 } );
	transform( ins, std::back_inserter( jumps ), []( const Audio * a ){ return a->getLength(); } );

	// Sum jumps to get mix positions
	std::vector<Second> startTimes;
	std::partial_sum( jumps.begin(), jumps.end() - 1, std::back_inserter( startTimes ), [&]( Second a, Second b ){ 
		return a + b + offset; 
		} );

	return mix( ins, startTimes, {} );
	}

Audio Audio::join( const std::vector<Audio> & ins, Second offset )
	{
	return join( getPtrs( ins ), offset );
	}

Audio Audio::select( const std::vector<const Audio *> & ins, const Func1x1 & selection, const std::vector<Second> & startTimes )
	{
	flan_FUNCTION_LOG;

	// Generate balances from selection
	std::vector<Function<Second, Amplitude>> balances;
	for( Index i = 0; i < ins.size(); ++i )
		{
		balances.emplace_back( [&selection, i]( Second t )
			{
			const float distance = std::abs( selection( t ) - i );
			if( distance >= 1 ) return 0.0f;
			else return std::sqrt( 1.0f - distance );
			} );
		}
		
	return mix( ins, startTimes, getPtrs( balances ) );
	}

Audio Audio::select( const std::vector<Audio> & ins, const Func1x1 & selection, std::vector<Second> startTimes )
	{
	return select( getPtrs( ins ), selection, startTimes );
	}