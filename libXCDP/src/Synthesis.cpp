#include "Synthesis.h"

#include <iostream>

#include <algorithm>

#include "WDL/resample.h"

static const double pi = std::acos( -1 );

namespace xcdp {
namespace Synthesis {

Audio waveform( RealFunc wave, double length, RealFunc freq, size_t samplerate, size_t oversample )
	{
	length = std::max( length, 0.0 );

	Audio::Format format;
	format.numChannels = 1;
	format.numSamples = length * samplerate;
	format.sampleRate = samplerate;
	Audio out( format );
	WDL_Resampler rs;
	rs.SetMode(true, 1, true, 512); // CPU-heavy but pretty good sinc interpolation
	double overrate = samplerate * oversample;
	rs.SetRates( overrate, samplerate );
	WDL_ResampleSample * rsinbuf = nullptr;
	int wanted = rs.ResamplePrepare( format.numSamples, 1, &rsinbuf );
	double phase = 0.0;
	
	for( size_t sample = 0; sample < wanted; ++sample )
		{
		double s = wave(phase);
		rsinbuf[sample] = s * 0.9; // lower the gain a bit in case the resampling causes overshoots
		phase += freq( double(sample) / overrate ) * ( 2 * pi / overrate );
		phase = std::fmod( phase, 2 * pi );
		
		}
	rs.ResampleOut( out.getBuffer()->data(), wanted, format.numSamples, 1 );
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