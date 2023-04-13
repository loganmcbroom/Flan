#include "flan/Audio/Audio.h"

#include <random>
#include <ctime>
#include <iostream>
#include <execution>
#include <ranges>

#include "flan/defines.h"

using namespace flan;

Audio Audio::texture( Time length, const Func1x1 & eventsPerSecond, const Func1x1 & scatter, const AudioMod & mod, bool feedback ) const
	{
	flan_PROCESS_START( Audio() );

	// Get constants
	const float secondsPerFrame = frameToTime();
	const float framesPerSecond = timeToFrame();
	const Frame lengthFrames = length * timeToFrame();

	// Input validation
	if( length <= 0 ) return Audio();
	auto eventsPerSecondSamples = eventsPerSecond.sample( 0, length * timeToFrame(), frameToTime() );
	std::for_each( eventsPerSecondSamples.begin(), eventsPerSecondSamples.end(), []( float & EpS )
		{ 
		EpS = std::max( 0.0f, EpS );
		} );
	auto eventsPerFrame = [&]( Frame f )
		{
		return eventsPerSecondSamples[ f ] * secondsPerFrame;
		};

	auto scatterSamples = scatter.sample( 0, length * timeToFrame(), frameToTime() );
	std::for_each( scatterSamples.begin(), scatterSamples.end(), []( float & scatter )
		{ 
		scatter = std::max( 0.0f, scatter );
		} );


	// Generate unscattered events by integrating eventsPerFrame
	// When the integral passes an integer, an event occurs
	std::vector<Frame> eventFrames;
	float eventAccumulator = 1.0f;
	for( Frame frame = 0; frame < lengthFrames; ++frame )
		{
		eventAccumulator += eventsPerFrame( frame );
		if( eventAccumulator >= 1.0f )
			{
			eventFrames.push_back( frame );
			eventAccumulator -= std::floor( eventAccumulator );
			}
		}
	if( eventFrames.empty() ) return Audio();

	// Scatter event frames, removing those landing outside lengthFrames
	std::default_random_engine rng( std::time( nullptr ) );
	std::for_each( std::execution::seq, eventFrames.begin(), eventFrames.end(), [&]( Frame & eventFrame )
		{
		// Scatter is standard deviation in events.
		const Time scatter_t = scatterSamples[ eventFrame ];
		const float eventsPerSecond_t = eventsPerSecondSamples[ eventFrame ];
		if( scatter_t == 0 ) return; // If there is no scatter, there is nothing to do
		if( eventsPerSecond_t == 0 ) return; // If there is no EpS, there is nothing to do

		const Time SDInSeconds = scatter_t / eventsPerSecond_t;
		const float	SDInFrames = SDInSeconds * framesPerSecond;
		const Frame scatteredFrame = std::normal_distribution<Time>( eventFrame, SDInFrames )( rng );
		if( 0 <= scatteredFrame && scatteredFrame < lengthFrames )
			eventFrame = scatteredFrame;
		else eventFrame = std::numeric_limits<Frame>::max(); // Using max to indicate "don't use", it will be sorted to the end shortly
		} );

	std::sort( std::execution::par_unseq, eventFrames.begin(), eventFrames.end() ); // Get our ducks back in a row

	// Shift event frames to align first event with 0 frame
	const Frame firstEventFrame = eventFrames[0];
	std::for_each( std::execution::par_unseq, eventFrames.begin(), eventFrames.end(), 
		[firstEventFrame]( Frame & f ){ f -= firstEventFrame; } );

	// Convert event frames back to times
	std::vector<Time> eventTimes( eventFrames.size() );
	std::transform( std::execution::par_unseq, eventFrames.begin(), eventFrames.end(), eventTimes.begin(), 
		[secondsPerFrame]( Frame frame ){ return frame * secondsPerFrame; } );

	// We need to handle all mod/no mod, and feedback/no feedback cases
	// If there is no mod feedback has no effect
	if( mod )
		{	
		std::vector<Audio> moddedAudio;
		moddedAudio.reserve( eventTimes.size() );
		if( feedback )
			{
			moddedAudio.push_back( mod( *this, 0 ) );
			// Sequential execution policy is required by feedback
			std::transform( eventTimes.begin() + 1, eventTimes.end(), std::back_inserter( moddedAudio ), [&]( Time eventTime )
				{
				const Audio & input = moddedAudio.back();
				return mod( input, eventTime );
				} );
			}
		else // Convert this using mod at each event time according to the mod execution policy
			{
			moddedAudio.resize( eventTimes.size() );
			runtimeExecutionPolicyHandler( mod.getExecutionPolicy(), [&]( auto policy ) {
				std::transform( policy, eventTimes.begin(), eventTimes.end(), moddedAudio.begin(), [&]( const Time eventTime ) {
					return mod( *this, eventTime );
					} ); 
				} );
			}
		return mix( moddedAudio, eventTimes, {} );
		}
	else // No mod, feedback does not apply
		{
		return mix( std::vector<const Audio *>( eventTimes.size(), this ), eventTimes, {} );
		}	
	}

Audio Audio::textureSimple( float length, const Func1x1 & eventsPerSecond, const Func1x1 & scatter, 
	const Func1x1 & repitch, const Func1x1 & gain, const Func1x1 & cutEnd, const Func1x1 & pan ) const
	{
	return texture( length, eventsPerSecond, scatter, [&]( const Audio & in, Time t )
		{
		Audio out = in.cut( 0, cutEnd( t ) )
			.repitch( repitch( t ) );
		out.modifyVolumeInPlace( gain( t ) );
		out.panInPlace( pan( t ) );
		return out;
		}, 
		false );
	}