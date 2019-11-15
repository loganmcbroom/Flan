#include <iostream>

#include <Windows.h>
#include <Mmsystem.h>

#include "Audio.h"
#include "PVOC.h"
#include "RealFunc.h"

using namespace xcdp;

void play( const Audio & toPlay )
	{
	toPlay.save( "TempFileSave.wav" );
	std::cout << "Playing sound ... ";
	if( ! PlaySound("TempFileSave.wav", nullptr, SND_FILENAME) )
		std::cout << "Error playing sound\n";
	std::cout << "Done\n";
	}

void main()
	{
	Audio bah = Audio( "Audio/bah.wav" )
	//.cut( 0.0, 0.2 )
	//.delay( .5, 10, .5 )
	//.fades( .5 )
	.setVolume( 0.8 );

	play( bah );

	//auto f = ADSR( .5, 2, 1, 3, .5, 
	//	.5, 2, 2 );
	//f.graph( "temp.tga" );
	//system( "start \"C:\\Program Files (x86)\\FastStone Image Viewer\\FSViewer.exe\" temp.tga " );
	}