#include "flan/Audio/Audio.h"

#include <iostream>

#include <algorithm>

#include "r8brain/CDSPResampler.h"

namespace flan {

Audio Audio::synthesize( const Func1x1 & wave, Time length, const Func1x1 & freq, size_t samplerate, size_t oversample )
	{
	flan_FUNCTION_LOG;

	if( length <= 0 ) return Audio();

	// Set up output
	Audio::Format format;
	format.numChannels = 1;
	format.numFrames = length * samplerate;
	format.sampleRate = samplerate;
	Audio out( format );

	const Frame numInFrames = out.getNumFrames() * oversample;
	const double inSamplerate = samplerate * oversample;

	auto freqSamples = freq.sample( 0, numInFrames, 1.0f / inSamplerate );

	// Get the phases that waveSamples needs to be evaluated at
	std::vector<float> phases( freqSamples.size() );
	std::exclusive_scan( std::execution::par_unseq, freqSamples.begin(), freqSamples.end(), phases.begin(), 0.0f, [&]( float a, float b ) -> float { 
		return std::fmod( a + b, inSamplerate ); } ); // Modular sum is associative, prefer exclusive_scan over partial_sum
	std::for_each( std::execution::par_unseq, phases.begin(), phases.end(), [&]( float & phase ){ phase /= inSamplerate; } ); // freq sum -> phase

	// Evaluate wave at phases using wave execution policy
	std::vector<float> waveSamples( phases.size() );
	runtimeExecutionPolicyHandler( wave.getExecutionPolicy(), [&]( auto policy ){
		std::transform( policy, phases.begin(), phases.end(), waveSamples.begin(), [&]( const float phase ){ return wave( phase ); } );
		} );

	// Downsample to requested sample rate
	r8b::CDSPResampler resampler( inSamplerate, out.getSampleRate(), waveSamples.size() );
	resampler.oneshot<float, float>( &waveSamples[0], waveSamples.size(), &out.getBuffer()[0], out.getBuffer().size() );

	return out;
	}

}
