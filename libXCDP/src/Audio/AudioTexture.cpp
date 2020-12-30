#include "xcdp/Audio.h"

#include <random>
#include <ctime>
#include <iostream>

#include "xcdp/defines.h"

using namespace xcdp;

// Concept for composable Audio Mods
//class AudioMod
//{
//public:
//	template< class T >
//	AudioMod( T )
//		: f( T )
//		{
//		}
//
//	Audio operator()( const Audio & in, float t ) const
//		{
//		return f( in, t );
//		}
//
//	AudioMod operator()( AudioMod other ) const
//		{
//		return AudioMod( [this, other]( const Audio & in, float t )
//			{
//			return operator()( other( in, t ), t );
//			} );
//		}
//
//	AudioMod operator||( AudioMod other ) const
//		{
//		return operator()( other );
//		}
//
//private:
//	std::function< Audio ( const Audio &, float ) > f;
//};



Audio Audio::texture( Time length, Func1x1 eventsPerSecond, Func1x1 scatter, Mod mod,
	bool feedback, XCDP_CANCEL_ARG_CPP ) const
	{
	XCDP_PROCESS_START( Audio() );

	// Input validation
	if( length <= 0 ) return Audio();
	Func1x1 safeEpS = [&eventsPerSecond]( Time t )
		{
		return std::max( 0.0f, eventsPerSecond( t ) );
		};
	Func1x1 safeScatter = [&scatter]( Time t )
		{
		return std::max( 0.0f, scatter( t ) );
		};

	// Get constants
	const float secondsPerFrame = frameToTime();
	const float framesPerSecond = timeToFrame();
	const Frame lengthFrames = length * timeToFrame();
	auto safeEpF = [&safeEpS, secondsPerFrame, framesPerSecond]( Frame f )
		{
		return safeEpS( f * secondsPerFrame ) * secondsPerFrame;
		};

	// Generate undistributed events by integrating eventsPerFrame
	// When the integral passes an integer, an event occurs
	std::vector<Frame> eventFrames;
	float eventAccumulator = 1.0f;
	for( Frame frame = 0; frame < lengthFrames; ++frame )
		{
		eventAccumulator += safeEpF( frame );
		if( eventAccumulator >= 1.0f )
			{
			XCDP_CANCEL_POINT( Audio() );
			eventFrames.push_back( frame );
			eventAccumulator -= std::floor( eventAccumulator );
			}
		}
	if( eventFrames.empty() ) return Audio();

	// Scatter event frames, removing those landing outside lengthFrames
	std::vector<Frame> scatteredEventFrames;
	std::default_random_engine rng( std::time( nullptr ) );
	for( auto & eventFrame : eventFrames )
		{
		// Scatter is standard deviation in events.
		const Time t = eventFrame * frameToTime();
		if( safeScatter( t ) == 0 )
			{
			scatteredEventFrames.push_back( eventFrame );
			continue;
			}
		if( safeEpS( t ) == 0 ) continue;
		const Time SDInSeconds = safeScatter( t ) / safeEpS( t );
		const float	SDInFrames = SDInSeconds * framesPerSecond;
		const Frame scatteredFrame = std::normal_distribution<Time>( eventFrame, SDInFrames )( rng );
		if( 0 <= scatteredFrame && scatteredFrame < lengthFrames )
			scatteredEventFrames.push_back( scatteredFrame );
		}
	if( scatteredEventFrames.empty() ) return Audio();
	std::sort( scatteredEventFrames.begin(), scatteredEventFrames.end() );

	// Shift event frames to align first event with 0 frame
	const Frame firstEventFrame = scatteredEventFrames[0];
	for( auto & event : scatteredEventFrames )
		event -= firstEventFrame;

	// Convert event frames back to times
	std::vector<Time> eventTimes( scatteredEventFrames.size() );
	std::transform( scatteredEventFrames.begin(), scatteredEventFrames.end(), eventTimes.begin(), 
		[secondsPerFrame]( Frame frame ){ return frame * secondsPerFrame; } );

	// Get mod outputs
	std::vector<Audio> moddedAudio;
	moddedAudio.reserve( eventTimes.size() );
	moddedAudio.push_back( mod? mod( *this, 0 ) : *this );
	for( int event = 1; event < eventTimes.size(); ++event )
		{
		XCDP_CANCEL_POINT( Audio() );
		const Time eventTime = eventTimes[event];
		const Audio & input = feedback ? moddedAudio[event - 1] : *this;
		moddedAudio.push_back( mod ? mod( input, eventTime ) : input );
		}

	return mix( moddedAudio, {}, eventTimes, canceller );
	}

Audio Audio::textureSimple( float length, Func1x1 eventsPerSecond, Func1x1 scatter, 
	Func1x1 repitch, Func1x1 gain, Func1x1 cutEnd, Func1x1 pan, XCDP_CANCEL_ARG_CPP ) const
	{
	return texture( length, eventsPerSecond, scatter, [&repitch, &gain, &cutEnd, &pan]( const Audio & in, Time t )
		{
		return in
			.cut( 0, cutEnd( t ) )
			.repitch( repitch( t ) )
			.modifyVolume( gain( t ) )
			.pan( pan( t ) );
		}, 
		false, canceller );
	}