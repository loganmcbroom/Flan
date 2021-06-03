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

Todo:
	pvoc generate from harmonics
		takes base freq as func, length, and a harmonic function
Task:
*/

#include <iostream>
#include <algorithm>
#include <numeric>

#include "xcdp/Audio.h"
#include "xcdp/PVOC.h" 
#include "xcdp/Function.h"
#include "xcdp/Wavetable.h"
#include "xcdp/Synthesis.h"
#include "xcdp/Graph.h"
#include "xcdp/FFTHelper.h"
#include "xcdp/WindowFunctions.h"
#include "xcdp/Utility/Timer.h"

using namespace xcdp;

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
	auto sine = Synthesis::sine( 1, []( float t ){ return 220 + 220 * t; } );
	auto saw = Synthesis::saw( .5, []( float t ){ return 220; } );

	auto p = yeahh.convertToPVOC();

	auto b = p.prism( []( Time t, int i, Frequency f, const std::vector<Magnitude> & ms ) 
		{ 
		return PVOC::MF{ ms[i], f * (i+1) * ( 1 + atan( t*3 ) * 0.1f * sin( t * 80 ) ) }; 
		}, true );

	graph( b.convertToGraph() );
	play( b.convertToAudio().setVolume( .7 ) );
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