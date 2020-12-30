#include "xcdp/Synthesis.h"

#include <iostream>

#include <algorithm>

#include "WDL/resample.h"

namespace xcdp {

static const float pi = std::acos( -1.0f );
static const float pi2 = 2.0f * pi;

namespace Synthesis {

Audio waveform( Func1x1 wave, Time length, Func1x1 freq, size_t samplerate, size_t oversample )
	{
	if( length <= 0 ) return Audio();

	// Set up output
	Audio::Format format;
	format.numChannels = 1;
	format.numFrames = length * samplerate;
	format.sampleRate = samplerate;
	Audio out( format );

	// Set up resampler
	WDL_Resampler rs;
	rs.SetMode( true, 1, true, 512 ); // CPU-heavy but pretty good sinc interpolation
	const double overrate = samplerate * oversample;
	rs.SetRates( overrate, double( samplerate ) );
	WDL_ResampleSample * rsinbuf = nullptr;
	const int wanted = rs.ResamplePrepare( format.numFrames, 1, &rsinbuf );

	float phase = 0.0f;
	for( Frame frame = 0; frame < wanted; ++frame )
		{
		const float s = wave( phase );
		rsinbuf[frame] = s * 0.9; // Lower the gain a bit in case the resampling causes overshoots
		phase += freq( float(frame) / overrate ) / overrate;
		phase = std::fmod( phase, 1.0f );
		}

	rs.ResampleOut( out.getSamplePointer( 0, 0 ), wanted, format.numFrames, 1 );

	return out;
	}

Audio sine( Time length, Func1x1 freq )
	{
	return waveform( []( float t ){ return std::sin( pi2 * t ); }, length, freq );
	}

Audio square( Time length, Func1x1 freq )
	{
	return waveform( []( float t ){ return t < 0.5f ? -1.0f : 1.0f; }, length, freq );
	}

Audio saw( Time length, Func1x1 freq )
	{
	return waveform( []( float t ){ return -1.0f + 2.0f * t; }, length, freq );
	}

Audio triangle( Time length, Func1x1 freq )
	{
	return waveform( []( float t ){ return t < 0.5f ? -1.0f + 4.0f * t : 3.0f - 4.0f * t; }, length, freq );
	}
}
}