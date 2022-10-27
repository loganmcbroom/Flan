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

int main()
	{
	auto bah = Audio( "Audio/Bah.wav" );
	auto bai = Audio( "Audio/Bai.wav" );
	auto meow = Audio( "Audio/meow.wav" );
	auto slaw = Audio( "Audio/slaw.wav" );
	auto sine = Synthesis::sine( 5, []( float t ){ return 220 + 220 * t; } );
	auto saw = Synthesis::saw( 5, []( float t ){ return 220 + 220 * t; } );

	auto b = bah.convertToPVOC().resonate( 0, .5 ).convertToAudio();


	play( b );
	//graph( b.convertToGraph() );

	return 0;
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