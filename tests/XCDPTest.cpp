

/* TODO:

Compatibility:
	check saving and loading on a big endian system

Additions:
	wavetable synth
	additive synth
	more filters
	add Func1x1 spline generator

Improvements:
	see if the Funct2x1 constructor can get working
	optimize cpu audio/pvoc conversions
	sort PVOC data by frequency, see if output improves
	stop documentation on wdl

Fixes:
	Make procs uncrashable
	test all procs after so many changes
	Switch to signed integer representation in PVOC RIFF data
	add interpolation to gpu methods with interpolator sampling
	stretch data loss? 0 initial frequency artifact?
	Accept saving 0 size audio files
	allow big widths in Audio::graph
	8bit support for loading

Dependancy removal:
	Add other lib compilation flags (implement)
	Add header only fft to build as option when fftw isn't used
	Write custom wav loader so libsndfile isn't required

check if the #ifdef in pvoc header can be removed
rewrite PVOCBuffer::save with new RIFF writer
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
	//Audio one = Synthesis::sine( 1, 440 );
	//one.save( "one.wav" );

	Audio two( "8bitSine.wav" );
	play( two.setVolume( .2 ) );
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