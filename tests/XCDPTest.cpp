/* TODO:

Compatibility:
	check saving and loading on a big endian system

Additions:
	wavetable synth
	additive synth
	more filters
	add Func1x1 spline generator

Improvements:
	see if the Funct2x1 fancy constructor can get working
	optimize cpu audio/pvoc conversions
	sort PVOC data by frequency, see if output improves
	stop documentation on wdl, etc

Fixes:
	Make procs uncrashable
	test all procs after so many changes
	add interpolation to gpu methods with interpolator sampling
	stretch data loss? 0 initial frequency artifact?
	allow big widths in Audio::graph
	check if \cond removed convertToAudio_FFTHelper doc
	move spline and djfft into some other folder?
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
	Audio audio = Synthesis::sine( .3, 440 ).setVolume( .9 );
	PVOC one = audio.convertToPVOC();
	one.save( "one.pvoc" );
	PVOC two( "one.pvoc" );

	std::cout << one << two;

	play( audio );
	play( two.convertToAudio().setVolume( .9 ) );
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