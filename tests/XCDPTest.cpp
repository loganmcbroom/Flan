

/* TODO:
wavetable synth
additive synth
more filters

stretch data loss? 0 initial frequency artifact?
sort PVOC data by frequency, see if output improves
Accept saving 0 size audio files
writebmp shouldn't take vector of vector (and fix documentation)
optimize cpu audio/pvoc conversions
add interpolation to gpu methods with interpolator sampling
Add more interpolators
Do you need to explicitly write out the forwarding constructors in Audio/PVOC?
Write custom wav loader so libsndfile isn't required
test all procs after so many changes
Add other lib compilation flags?
Make procs uncrashable
Implement Audio::fade flags
Move timer to utilities
Switch to signed integer representation in PVOC RIFF data
see if the Funct2x1 constructor can get working
add Func1x1 spline generator?
*/

#include <iostream>

#include <algorithm>

#include "xcdp/Audio.h"
#include "xcdp/PVOC.h"
#include "xcdp/Function.h"
#include "xcdp/Synthesis.h"

using namespace xcdp;

void play( const Audio & toPlay );
void graph( const std::string & f );

void main()
	{

	Func1x1 func = []( float t ){ return t*t; };
	func.periodize( []( float t ){ return 1.0f/(t*t+1.0);} )
	.graph( "tempgraph.bmp" );

    //Audio in = Audio( "Audio/meow.wav" );

	//auto out = 
    //in
	//.convertToPVOC()
 //   .modify( []( vec2 v ) 
 //       { 
 //       v.y() += ( 1.0 - cos( v.x() * 10.0f ) ) * 100.0f;
 //       v.x() += ( 1.0f - cos( v.x() * 10.0f ) );
 //       return v;
 //       })
 //   .graph( "tempgraph.bmp" )
	//.convertToAudio()
 //   ;

    graph( "tempgraph.bmp" );
	//play( out.setVolume( .5 ) );

 //   PVOC meow = Audio( "Audio/meow.wav" ).convertToPVOC();

 //   Timer timer;
 //   timer.start();
	//for( int i = 0; i < 100; ++i )
	//	 meow.modifyFrequency( []( double t, double f ){ return t; } );
 //   timer.stop();
 //   std::cout << timer.elapsedSeconds() << std::endl;
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