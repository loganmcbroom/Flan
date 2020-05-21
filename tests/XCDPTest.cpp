#include <iostream>

#include <algorithm>

#include "xcdp/Audio.h"
#include "xcdp/PVOC.h"
#include "xcdp/Function.h"
#include "xcdp/Synthesis.h"

using namespace xcdp;

void play( const Audio & toPlay );
void graph( const std::string & f );

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
test gitignore inside a folder so it doesn't get everything if it doesn't work
test all procs after so many changes
Add WDL compilation flag
Add other compilation flags?
Make procs uncrashable
Implement Audio::fade flags
Move timer to utilities
Switch to signed integer representation in PVOC RIFF data
see if the Funct2x1 constructor can get working
add Func1x1 spline generator?
fix oscillate, it doesn't actually loop, + doc
*/

#include <string>
#include <windows.h>

std::string getexepath()
{
  char result[ MAX_PATH ];
  return std::string( result, GetModuleFileName( NULL, result, MAX_PATH ) );
}

class Timer
{
public:
    void start()
    {
        m_StartTime = std::chrono::system_clock::now();
        m_bRunning = true;
    }
    
    void stop()
    {
        m_EndTime = std::chrono::system_clock::now();
        m_bRunning = false;
    }
    
    float elapsedMilliseconds()
    {
        std::chrono::time_point<std::chrono::system_clock> endTime;
        
        if(m_bRunning)
        {
            endTime = std::chrono::system_clock::now();
        }
        else
        {
            endTime = m_EndTime;
        }
        
        return std::chrono::duration_cast<std::chrono::milliseconds>(endTime - m_StartTime).count();
    }
    
    float elapsedSeconds()
    {
        return elapsedMilliseconds() / 1000.0;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> m_StartTime;
    std::chrono::time_point<std::chrono::system_clock> m_EndTime;
    bool                                               m_bRunning = false;
};

void main()
	{

    Audio in = Audio( "Audio/meow.wav" );

	auto out = 
    in
	.convertToPVOC()
    .modify( []( vec2 v ) 
        { 
        v.y() += ( 1.0 - cos( v.x() * 10.0f ) ) * 100.0f;
        v.x() += ( 1.0f - cos( v.x() * 10.0f ) );
        return v;
        })
    .graph( "tempgraph.bmp" )
	.convertToAudio()
    ;

    graph( "tempgraph.bmp" );
	play( out.setVolume( .5 ) );

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