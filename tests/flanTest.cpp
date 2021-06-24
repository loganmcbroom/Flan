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

	auto p = bah.convertToPVOC();

	//auto b = PVOC::generate( 1, 111, []( float t, float i ){ return 1 / ( 1 + i ); } );

	//auto b = p.prism( []( Time t, int i, Frequency f, const std::vector<Magnitude> & ms ) 
	//	{ 
	//	return PVOC::MF{ ms[i], (75.0f + 40 * t ) * (i+1) }; 
	//	}, false );

	//play( b.convertToAudio() );
	//play( irr );
	//graph( irr.convertToPVOC().convertToGraph() );
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