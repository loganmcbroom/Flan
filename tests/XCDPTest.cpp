#include <iostream>

#include <algorithm>

#include "Audio.h"
#include "PVOC.h"
#include "RealFunc.h"

using namespace xcdp;

void play( const Audio & toPlay );

//TODO
//solve stretch missing frame issue
//Deal with interpolate upper bound nonsense
//find a valid convolve?
//finish accumulate
//sort PVOC data by frequency, see if output improves
//qt

Audio makeSine( RealFunc freq )
	{
	Audio::Format format;
	format.numChannels = 1;
	format.numSamples = 48000;
	format.sampleRate = 48000;
	Audio out( format );

	for( size_t sample = 0; sample < out.getNumSamples(); ++sample )
		out.setSample( 0, sample, std::sin( freq( out.sampleToTime( sample ) ) * 3.14159 * 2.0 * double( sample ) / double( out.getSampleRate() ) ) );
	
	return out;
	}

void main()
	{
	//Audio sine = makeSine( [](double t){ return 440.0*(1+t);} );
	//PVOC meow = Audio( "Audio/meow.wav" ).convertToPVOC();
	Audio bah = Audio( "Audio/Beep.wav" )
	.convertToPVOC( 1024, 8 )
	.interpolate( .03 )
	.convertToAudio()
	.setVolume( 0.8 );

	RealFunc tvp1 = RealFunc::interpolatePoints( { {0, 2}, {1, 2}, {2, 1}, {3, 2} } );
	RealFunc tvp2 = RealFunc::interpolatePoints( { {0, 2}, {1, 2}, {2, 1}, {3, 2} }, Interpolators::sine );
	//tvp1.graph();

	//system( "start \"C:\\Program Files (x86)\\FastStone Image Viewer\\FSViewer.exe\" tempRealFuncGraph.tga" );
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