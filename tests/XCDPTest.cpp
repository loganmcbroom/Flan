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
	std::cout << "Playing sound ... \n";
	if( ! PlaySound("TempFileSave.wav", nullptr, SND_FILENAME) )
		std::cout << "Error playing sound\n";
	}

void main()
	{
	PVOC meowPVOC = Audio( "Audio/meow.wav" ).convertToPVOC();
	Audio bah = Audio( "Audio/bah.wav" )
	.convertToPVOC()
	.replaceAmplitudes( meowPVOC )
	.getSpectrograph( "temp.tga" )
	.convertToAudio()
	.setVolume( 0.8 );


	//auto f = ADSR( .5, 2, 1, 3, .5, 
	//	.5, 2, 2 );
	//f.graph( "temp.tga" );

	system( "start \"C:\\Program Files (x86)\\FastStone Image Viewer\\FSViewer.exe\" temp.tga " );
	play( bah );
	}