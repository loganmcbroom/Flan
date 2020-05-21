#include "xcdp/Synthesis.h"

#include <iostream>

#include <algorithm>

#include "WDL/resample.h"

static const float pi = std::acos( -1.0f );

namespace xcdp {
namespace Synthesis {

Audio waveform( Func1x1 wave, float length, Func1x1 freq, size_t samplerate, size_t oversample )
	{
	length = std::max( length, 0.0f );

	Audio::Format format;
	format.numChannels = 1;
	format.numFrames = length * samplerate;
	format.sampleRate = samplerate;
	Audio out( format );
	WDL_Resampler rs;
	rs.SetMode(true, 1, true, 512); // CPU-heavy but pretty good sinc interpolation
	double overrate = samplerate * oversample;
	rs.SetRates( overrate, double( samplerate ) );
	WDL_ResampleSample * rsinbuf = nullptr;
	int wanted = rs.ResamplePrepare( format.numFrames, 1, &rsinbuf );
	float phase = 0.0;
	
	for( size_t sample = 0; sample < wanted; ++sample )
		{
		float s = wave(phase);
		rsinbuf[sample] = s * 0.9; // lower the gain a bit in case the resampling causes overshoots
		phase += freq( float(sample) / overrate ) * ( 2.0f * pi / overrate );
		phase = std::fmod( phase, 2 * pi );
		
		}
	rs.ResampleOut( out.getSamplePointer( 0, 0 ), wanted, format.numFrames, 1 );
	return out;
	}

Audio sine( float length, Func1x1 freq )
	{
	return waveform( []( float t ){ return std::sin( t ); }, length, freq );
	}
Audio square( float length, Func1x1 freq )
	{
	return waveform( []( float t ){ return t < pi ? -1.0 : 1.0; }, length, freq );
	}
Audio saw( float length, Func1x1 freq )
	{
	return waveform( []( float t ){ return -1.0 + 1.0 / pi * t; }, length, freq );
	}
Audio triangle( float length, Func1x1 freq )
	{
	return waveform( []( float t ){ return t < pi? -1.0 + 2.0 / pi * t : 3.0 - 2.0 / pi * t; }, length, freq );
	}
}
}