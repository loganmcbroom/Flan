#include <iostream>

#include <algorithm>

#include "Audio.h"
#include "PVOC.h"
#include "RealFunc.h"
#include "Synthesis.h"

using namespace xcdp;

void play( const Audio & toPlay );
void graph( const std::string & f );

//TODO
//AA filter synthesis
//desample ending bug
//sort PVOC data by frequency, see if output improves
//Make sure surfaces are supported everywhere they could be
//qt

void main()
	{
	//Audio sine = makeSine( [](double t){ return 440.0*(1+t);} );
	Audio other = Audio( "Audio/bah.wav" );
	Audio bah = Audio( "Audio/co.wav" )
	//.setVolume( 0.5 )
	//Audio bah = Synthesis::saw( 2.0, []( double t ){ return 220.0+220.0*t; } )
	//.convertToPVOC( 2048, 32 )
	//.modifyFrequency( []( double t, double f ){ return f*(0.7 + 0.5*sin(16.0*t)); } )
	//.modifyTime( []( double t, double f ){ return t*1.0 + f/12000.0; } )
	//.resonate( 8, .3 )
	//.getSpectrograph( "temp.bmp" )
	//.convertToAudio()
	.convolve( other )
	.setVolume( 0.7 );

	//RealFunc tvp1 = RealFunc::interpolatePoints( { {0, 2}, {1, 2}, {2, 1}, {3, 2} } );
	//RealFunc tvp2 = RealFunc::interpolatePoints( { {0, 2}, {1, 2}, {2, 1}, {3, 2} }, Interpolators::sine );
	//tvp1.graph( "tempBMP.bmp" );
	
	graph( "temp.bmp" );
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