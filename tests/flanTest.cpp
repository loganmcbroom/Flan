/* Ideas:
	filter package - just do it in PVOC?
	continuity maintenance in wavelength/frequency finders
	octave error minimization in yin-based pitch finder (fast contour finder)
	audio to function of frequency
	convolve constant signals via fft
	Add \mainpage with examples somewhere
	Audio::removeSilence
	Audio::cut uses an extra buffer copy when calling fades, remove it
	contour pitch mean should be magnitude weighted
	only generate docs for flan stuff
	pvoc::alignHarmonics

Todo:

Task:
	select seems slow, check it out
	select/mix/join should consider bitrate conversion when needed
	delete and rebuild sndfile to make sure it's in release mode

	add "list( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/")" to cmake.in
*/

#include <iostream>
#include <algorithm>
#include <numeric>

#include "flan/Audio.h"
#include "flan/PVOC.h" 
#include "flan/Function.h"
#include "flan/Wavetable.h"
#include "flan/Synthesis.h"
#include "flan/Graph.h"
#include "flan/FFTHelper.h"
#include "flan/WindowFunctions.h"
#include "flan/Utility/Timer.h"

using namespace flan;

void play( const Audio & toPlay );
void graph( const std::string & f );
void graph( const bitmap_image & b );

void main()
	{
	auto bah = Audio( "Audio/Bah.wav" );
	auto bai = Audio( "Audio/Bai.wav" );
	auto meow = Audio( "Audio/meow.wav" );
	auto slaw = Audio( "Audio/slaw.wav" );
	auto yeahh = Audio( "Audio/feel.wav" );
	auto sine = Synthesis::sine( 5, []( float t ){ return 220 + 220 * t; } );
	auto saw = Synthesis::saw( 5, []( float t ){ return 220 + 220 * t; } );

	auto b = bah.convertToPVOC().resonate( 0, .5 ).convertToAudio();


	play( b );
	//graph( b.convertToGraph() );
	}


#include <Windows.h>
#include <Mmsystem.h>
void play( const Audio & toPlay )
	{
	if( !toPlay.save( "TempFileSave.wav" ) )
        return;
	std::cout << "Playing sound ... \n";
	if( ! PlaySound("TempFileSave.wav", nullptr, SND_FILENAME) )
		std::cout << "Error playing sound\n";
	}

void graph( const std::string & f )
	{
	system( (std::string("start ") + f).c_str() );
	}

void graph( const bitmap_image & b )
	{
	b.save_image( "Temp.bmp" );
	system( (std::string("start ") + "Temp.bmp").c_str() );
	}