#include <iostream>
#include <algorithm>
#include <numeric>
#include <ranges>

#include "flan/Audio/Audio.h"
#include "flan/PV/PV.h" 
#include "flan/SPV/SPV.h" 
#include "flan/SQPV/SQPV.h"
#include "flan/Function.h"
#include "flan/Wavetable.h"
#include "flan/Synthesis.h"
#include "flan/Graph.h"
#include "flan/WindowFunctions.h"
#include "flan/Utility/Timer.h"
#include "Iir1/Iir.h"
#include "flan/DSPUtility.h"

using namespace flan;

void play( const Audio & toPlay );
void graph( const std::string & f );
void graph( const bitmap_image & b );
template<typename T, typename Callable>
void test( const T & a, int n, Callable c );

const float pi = std::acos( -1.0f );

int main()
	{
	auto synth = Audio::synthesize( []( float t ){ return t / flan::pi - 1.0f; }, 1, []( float t ){ return 200.0f; } );
	auto bah = Audio( "Bah.wav" );

	//play( bah.convertToPV().convertToAudio().setVolume( .7 ) );
	//play( bah.convertToSPV().convertToAudio().setVolume( .7 ) );
	play( bah.convertToMidSideSQPV().modifyFrequency( []( vec2 tf ){ return tf.f() * .5f; } ).convertToLeftRightAudio().setVolume( .7 ) );
	
	return 0;
	}


#include <Windows.h>
#include <Mmsystem.h>
void play( const Audio & toPlay )
	{
	if( !toPlay.save( "TempFileSave.wav" ) )
        return;
	std::cout << "Playing sound ... " << std::endl;
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

template<typename T, typename Callable>
void test( const T & a, int n, Callable c )
	{
	Timer t;
	t.start();
	for( const int i : std::views::iota( 0, n ) )
		c( a );
	t.stop();
	std::cout << t.elapsedMilliseconds() / n << " ms per call" << std::endl;
	}

