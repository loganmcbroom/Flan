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
analysis file saving
documentation
more filters
qt
*/

void main()
	{
	Audio( "Audio/bah.wav" )
	.convertToPVOC()
	.save( "spectral.pvoc" );
	
	Audio bah = PVOC( "spectral.pvoc" )
	.convertToAudio();

	
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
	const std::string command = std::string("start \"C:\\ProgramData\\Microsoft\\Windows\\Start Menu\\Programs\\paint.exe\" ");
	system( (command + f).c_str() );
	}