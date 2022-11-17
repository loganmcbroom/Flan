#include "flan/Audio.h"

#include <iostream>

#include <algorithm>

#include "WDL/resample.h"

namespace flan {

Audio Audio::synthesize( Func1x1 wave, Time length, Func1x1 freq, size_t samplerate, size_t oversample, flan_CANCEL_ARG_CPP )
	{
	flan_FUNCTION_LOG;

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
		flan_CANCEL_POINT( Audio() );
		const float s = wave( phase );
		rsinbuf[frame] = s * 0.9; // Lower the gain a bit in case the resampling causes overshoots
		phase += freq( float(frame) / overrate ) / overrate;
		phase = std::fmod( phase, 1.0f );
		}

	rs.ResampleOut( out.getSamplePointer( 0, 0 ), wanted, format.numFrames, 1 );

	return out;
	}

//Audio waveform( Func1x1 wave, Time length, Func1x1 freq, size_t samplerate, size_t oversample )
//	{
//	if( length <= 0 ) return Audio();
//
//	// Set up output
//	Audio::Format format;
//	format.numChannels = 1;
//	format.numFrames = length * samplerate;
//	format.sampleRate = samplerate;
//	Audio out( format );
//
//	// Set up resampler
//	WDL_Resampler rs;
//	rs.SetMode( true, 1, true, 512 ); // CPU-heavy but pretty good sinc interpolation
//	const double overrate = samplerate * oversample;
//	rs.SetRates( overrate, double( samplerate ) );
//	WDL_ResampleSample * rsinbuf = nullptr;
//	const int wanted = rs.ResamplePrepare( format.numFrames, 1, &rsinbuf );
//
//	float phase = 0.0f;
//	for( Frame frame = 0; frame < wanted; ++frame )
//		{
//		const float s = wave( phase );
//		rsinbuf[frame] = s * 0.9; // Lower the gain a bit in case the resampling causes overshoots
//		phase += freq( float(frame) / overrate ) / overrate;
//		phase = std::fmod( phase, 1.0f );
//		}
//
//	rs.ResampleOut( out.getSamplePointer( 0, 0 ), wanted, format.numFrames, 1 );
//
//	return out;
//	}

}
