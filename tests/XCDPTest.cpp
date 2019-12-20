#include <iostream>

#include <algorithm>

#include "Audio.h"
#include "PVOC.h"
#include "RealFunc.h"
#include "Synthesis.h"

using namespace xcdp;

void play( const Audio & toPlay );
void graph( const std::string & f );

/* TODO:
stretch data loss? 0 initial frequency artifact?
sort PVOC data by frequency, see if output improves
Make sure surfaces are supported everywhere they could be

wavetable synth
additive synth
documentation
more filters
qt
*/

void main()
	{
	Audio bah = Audio( "Audio/Bah.wav" )
	.repitch( []( double t ){ return .3+t; } );
	
	//graph( "temp.bmp" );
	play( bah );
	}

#include <Windows.h>
#include <Mmsystem.h>
void play( const Audio & toPlay )
	{
	toPlay.save( "TempFileSave.wav" );
	std::cout << "Playing sound ... \n";
	if( ! PlaySound("TempFileSave.wav", nullptr, SND_FILENAME) )
		std::cout << "Error playing sound\n";
	}

void graph( const std::string & f )
	{
	system( (std::string("start ") + f).c_str() );
	}