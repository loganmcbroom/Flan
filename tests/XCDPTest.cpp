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
accept saving 0 size audio files

wavetable synth
additive synth
documentation
more filters
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

#include <complex>
void main()
	{
 //   Audio in = Audio( "Audio/meow.wav" );
	//Audio out = in
	//.convertToPVOC()
	//.convertToAudio();

	//play( out.setVolume( .5 ) );


    /*
    Audio -> PVOC
    16.791
    5.458 
    2.851

    PVOC -> Audio
    12.46
    7.562
    3.529
    */

    Audio meow( "Audio/meow.wav" );

    Timer timer;
    timer.start();
	for( int i = 0; i < 100; ++i )
		meow.convertToPVOC();
    timer.stop();
    std::cout << timer.elapsedSeconds() << std::endl;
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