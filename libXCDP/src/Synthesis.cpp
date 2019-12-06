#include "Synthesis.h"

#include <iostream>

#include <algorithm>

static const double pi = std::acos( -1 );

namespace xcdp {
namespace Synthesis {

Audio waveform( RealFunc wave, double length, RealFunc freq )
	{
	length = std::max( length, 0.0 );

	Audio::Format format;
	format.numChannels = 2;
	format.numSamples = length * 48000;
	format.sampleRate = 48000;
	Audio out( format );

	for( size_t channel = 0; channel < out.getNumChannels(); ++channel )
		for( size_t sample = 0; sample < out.getNumSamples(); ++sample )
			{
			const double input = 2.0 * pi * freq( out.sampleToTime( sample ) ) * double(sample) / out.getSampleRate();
			out.setSample( channel, sample, wave( std::fmod( input, 2 * pi ) ) );
			}

	
	
	return out;
	}

Audio sine( double length, RealFunc freq)
	{
	return waveform( [](double t){ return std::sin( t ); }, length, freq );
	}
Audio square( double length, RealFunc freq)
	{
	return waveform( [](double t){ return t < pi ? -1.0 : 1.0; }, length, freq );
	}
Audio saw( double length, RealFunc freq)
	{
	return waveform( [](double t){ return -1.0 + 1.0 / pi * t; }, length, freq );
	}
Audio triangle( double length, RealFunc freq)
	{
	return waveform( [](double t){ return t < pi? -1.0 + 2.0 / pi * t : 3.0 - 2.0 / pi * t; }, length, freq );
	}
}
}