/** TODO:
	housekeep chans
	extend loop
	filter userbank
	filter varibank(2)
	RealFunc ADSR( double A, double D, double S, double R, double Aexp = 0, double Dexp = 0, double Rexp = 0 );
	distort average
	distort interpolate
	distort multiply
*/

#include "Procs.h"

#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>

#include <samplerate.h>

namespace xcdp {

const double pi = std::acos( -1.0 );

//========================================================
// Utility Functions
//========================================================

// You want print to log? $5 print to log best deal in town.
void printToLog( std::string s )
	{
	std::cout << s << std::endl;
	}
bool doSampleRatesMatch( ProcInputVec ins )
	{
	//Check if all sample rates match the first file
	int sampleRate = ins[0].getSampleRate();
	for( auto & in : ins ) 
		if( in.getSampleRate() != sampleRate )
			{
			printToLog( "Mismatched sample rates" );
			return false;
			}
	return true;
	}
bool doChannelCountsMatch( ProcInputVec ins )
	{
	//Check if all channel counts match the first file
	int numChannels = ins[0].getNumChannels();
	for( auto & in : ins ) 
		if( in.getNumChannels() != numChannels )
			{
			printToLog( "Mismatched channel count" );
			return false;
			}
	return true;
	}

int getMaxNumChannels( ProcInputVec ins )
	{
	return std::max_element( ins.begin(), ins.end(), []( ProcInput a, ProcInput b )
		{ 
		return a.getNumChannels() < b.getNumChannels();
		} )->getNumChannels();
	}
int getMaxNumSamples( ProcInputVec ins )
	{
	return std::max_element( ins.begin(), ins.end(), []( ProcInput a, ProcInput b )
		{ 
		return a.getNumFrames() < b.getNumFrames();
		} )->getNumFrames();
	}

float getMaxSample( ProcInput in )
	{
	return std::abs( *std::max_element( 
			in.buffer.begin(), 
			in.buffer.end(), 
			[]( float a, float b ){ return std::abs(a) < std::abs(b); }));
	}

AudioBuffer mono2Stereo( ProcInput mono )
	{
	if( mono.getNumChannels() != 1 )
		printToLog( "That wasn't a mono file you tried to make stereo just now.\n" );

	AudioBuffer stereo;
	stereo.copyFormat( mono );
	stereo.setBufferSize( 2, mono.getNumFrames() );

	for( int frame = 0; frame < stereo.getNumFrames(); ++frame )
		{
		stereo.setSample( 0, frame, mono.getSample( 0, frame ) );
		stereo.setSample( 1, frame, mono.getSample( 0, frame ) );
		}

		return stereo;
	}

//========================================================
// Procs
//========================================================

AudioBuffer modifyVolume( ProcInput in, RealFunc volumeLevel )
	{
	AudioBuffer out;
	out.copyFormat( in );
	out.setBufferSize( in.getNumChannels(), in.getNumFrames() );

	for( int channel = 0; channel < in.getNumChannels(); ++channel )
		for( int frame = 0; frame < in.getNumFrames(); ++frame )
			{
			double calculatedSample = in.getSample( channel, frame )*volumeLevel(in.getTimeOfFrame(frame));
			out.setSample( channel, frame, calculatedSample );
			}

	return out;
	}
AudioBuffer modifyVolume( ProcInput in, double volumeLevel )
	{
	return modifyVolume( in, [ volumeLevel ]( double ){ return volumeLevel; } );
	}

AudioBuffer setVolume( ProcInput in, double level )
	{
	// Divide by maxVolume to normalize, multiply by level to set
	double volume = level / getMaxSample( in );
	return modifyVolume( in, [ volume ]( double t ){ return volume; } );
	}

AudioBuffer mix( ProcInputVec ins, std::vector< RealFunc > balances )
	{
	//If we didn't get a balances input, construct a default one, otherwise check input validity
	if( balances.size() == 0 )
		{
		double insSize = ins.size();
		auto defaultBalanceFunc = [insSize]( double )
			{ 
			return 1.0 / insSize; //1/n trades off natural sound for guaranteed no clip
			};
		balances.assign( ins.size(), defaultBalanceFunc );
		}
	else if( balances.size() < ins.size() )
		{
		printToLog( "Too few balances sent to mix, padding with 1/n's. \n" );
		double oneOverN = 1.0 / double( ins.size() ); 
		balances.insert( balances.end(), ins.size() - balances.size(), [oneOverN]( double )
			{ 
			return oneOverN;
			} );
		}
	else if( balances.size() > ins.size() )
		{
		printToLog( "Too Many balances sent to mix. Not really a problem, just thought you should know. \n" );
		}

	if( !doChannelCountsMatch( ins ) )
		{
		printToLog( "You can't mix files with differing channel counts at the moment (Is there a good way?).\n" );
		return AudioBuffer();
		}

	int numChannels = getMaxNumChannels( ins );
	int numSamples  = getMaxNumSamples ( ins );

	AudioBuffer out;
	out.copyFormat( ins[0] );
	out.setBufferSize( numChannels, numSamples );

	for( int channel = 0; channel < numChannels; ++channel )
		for( int frame = 0; frame < numSamples; ++frame )
			{
			double mixedSample = 0.0;
			for( int i = 0; i < ins.size(); ++i )
				mixedSample += balances[i]( ins[0].getTimeOfFrame( frame ) ) * ins[i].getSample( channel, frame );
			out.setSample( channel, frame, std::clamp( mixedSample, -1.0, 1.0 ) );
			}

	return out;
	}
AudioBuffer mix( ProcInputVec ins, const std::vector< double > & balances )
	{
	std::vector< RealFunc > funcs( balances.size() );
	for( unsigned int i = 0; i < balances.size(); ++i )
		funcs[i] = [i, &balances](double){ return balances[i]; };

	return mix( ins, funcs );
	}

AudioBuffer waveshape( ProcInput in, RealFunc shaper )
	{
	AudioBuffer out;
	out.copyFormat( in );
	out.setBufferSize( in.getNumChannels(), in.getNumFrames() );

	for( int channel = 0; channel < in.getNumChannels(); ++channel )
		for( int frame = 0; frame < in.getNumFrames(); ++frame )
			{
			out.setSample( channel, frame, shaper( in.getSample( channel, frame ) ) );
			}

	return out;
	}

AudioBuffer pan( ProcInput in, RealFunc panAmount )
	{
	//Stereo panning algorithm
	auto stereoPan = []( ProcInput in, RealFunc panAmount )
		{
		AudioBuffer out;
		out.copyFormat( in );
		out.setBufferSize( in.getNumChannels(), in.getNumFrames() );

		for( int channel = 0; channel < in.getNumChannels(); ++channel )
			for( int frame = 0; frame < in.getNumFrames(); ++frame )
				{
				double sample = in.getSample( channel, frame )
					* sin( pi / 4.0 * ( std::clamp( panAmount( in.getTimeOfFrame(frame)), -1.0, 1.0 ) + 3.0 - double(channel)*2.0 ) )
					* sqrt( 2.0 );
				out.setSample( channel, frame, sample );
				}
		return out;
		};
	
	switch( in.getNumChannels() )
		{
		case 1: //Panning mono: cast to stereo and do stereo pan
			return stereoPan( mono2Stereo( in ), panAmount );
			break;
		case 2: //Panning stereo: Use stereo pan as defined above
			return stereoPan( in, panAmount );
			break;
		default:
			printToLog( "I don't know how to pan that number of channels" );
			return in;
		}
	}
AudioBuffer pan( ProcInput in, double panAmount )
	{
	return pan( in, [panAmount](double){ return panAmount; } );
	}

AudioBuffer iterate( ProcInput in, int n, std::function< AudioBuffer (ProcInput, int n) > mod )
	{
	AudioBuffer out;
	out.copyFormat( in );
	out.setNumChannels( in.getNumChannels() );

	//For each iteration
	for( int i = 0; i < n; ++i )
		{
		//Either modify in with mod or if there was no mod, don't.
		//Should be fine speed wise with move copy
		const AudioBuffer inMod = mod == 0 ? in : mod( in, i );

		//Append inMod to out
		out.buffer.insert( out.buffer.end(), inMod.buffer.begin(), inMod.buffer.end() );
		out.info.frames += inMod.getNumFrames();
		}

	return out;
	}

AudioBuffer cutAtSamples( ProcInput in, int startSample, int endSample )
	{
	//input validity checking
	startSample = std::min( 0, startSample );
	endSample = std::max( in.getNumFrames(), endSample );
	if( endSample < startSample ) endSample = startSample;

	AudioBuffer out;
	out.copyFormat( in );
	out.setBufferSize( in.getNumChannels(), endSample - startSample );

	for( int frame = 0; frame < out.getNumFrames(); ++frame )
		for( int channel = 0; channel < out.getNumChannels(); ++channel )
			out.setSample( channel, frame, in.getSample( channel, startSample + frame ) );

	return out;
	}
AudioBuffer cutAtTimes( ProcInput in, double startTime, double endTime )
	{
	//Input time validity is checked in cutAtSamples
	return cutAtSamples( in, startTime / in.getSampleRate(), endTime / in.getSampleRate() );
	}

AudioBuffer join( AudioBufferVec ins )
	{
	int numChannels = getMaxNumChannels( ins );
	int numSamples = 0;
	for( auto & in : ins ) numSamples += in.getNumFrames();

	AudioBuffer out;
	out.copyFormat( ins[0] );
	out.setBufferSize( numChannels, numSamples );

	int currentOutStartSample = 0;

	//For each in, copy in into out
	for( int in = 0; in < ins.size(); ++in )
		{
		for( int frame = 0; frame < ins[in].getNumFrames(); ++frame )
			for( int channel = 0; channel < ins[in].getNumChannels(); ++channel )
				out.setSample( channel, currentOutStartSample + frame, ins[in].getSample( channel, frame ) );
		currentOutStartSample += ins[in].getNumFrames();
		}

	return out;
	}

AudioBuffer repitch( ProcInput in, RealFunc factor )
	{
	//We can't have true continuous pitch shifting, to we approximate the scaling factor
	// piecewise linearly. blocksize decides how often the continuous factor is sampled
	const int blockSize = in.getSampleRate() / 32; //32 times per second
	const int numBlocks = floor( in.getNumFrames() / blockSize );

	//Protection from bad factor outputs
	auto safeFactor = [factor]( double t ){ return std::clamp( factor(t), .00001, 100000.0 ); };

	//Estimate output buffer length via integration (trapezoidal approximation)
	//libsamplerate uses linear pitch shifting as well so this should be within a sample
	//of the actual needed buffer length. It is rounded up, and there is an additional 
	//protection against bad array access later.
	double outLength = 0.5 / safeFactor( 0 ) 
					 + 0.5 / safeFactor( in.getTimeOfFrame( blockSize * numBlocks ) ); // first and last terms
	for( int i = 1; i < numBlocks; ++i ) //Middle terms
		outLength += 1.0 / safeFactor( in.getTimeOfFrame( i * blockSize ) );
	outLength *= blockSize; //Multiply by trapezoid base size
	//Handle remaining ( in length modulo blockSize ) samples
	outLength += ( in.getNumFrames() % blockSize ) 
					* ( 0.5 / safeFactor( in.getTimeOfFrame( numBlocks * blockSize ) ) 
					  + 0.5 / safeFactor( in.getTimeOfFrame( in.getNumFrames() )   ) );

	AudioBuffer out;
	out.copyFormat( in );
	out.setBufferSize( in.getNumChannels(), ceil( outLength ) );

	//Converter initialization
	int err = 0;
	SRC_STATE * state = src_new( 1, in.getNumChannels(), &err ); //type 1 resampling
	if( state == nullptr )
		{
		printToLog( "Error creating sample rate converter in repitch.\n" );
		return in;
		}

	//Conversion settings
	SRC_DATA data;

	//What the heck are these two? These, and the block+1 and numBlocks+1 later, are 
	// there to correctly lerp between pitches. If we don't start at 0 and set the initial
	// pitch, the first block will have constant pitch
	data.input_frames  = 0; 
	src_process( state, &data );

	data.data_in = in.buffer.data(); //We will increment these as needed in the main loop
	data.data_out = out.buffer.data();
	data.input_frames  = blockSize; //Process one block at a time
	data.output_frames = out.getNumFrames();
	data.end_of_input = 0;

	//Process blocks
	for( int block = 0; block < numBlocks; ++block )
		{
		//Note the 1/factor, we are increasing pitch by decreasing sample rate and then playing bate at
		//the original sample rate
		data.src_ratio = 1.0 / safeFactor( in.getTimeOfFrame( (block + 1) * blockSize ) );
		src_process( state, &data );
		data.data_in  += blockSize * in.getNumChannels();
		data.data_out += data.output_frames_gen * out.getNumChannels();
		data.output_frames -= data.output_frames_gen; //This protects from bad array access
		}

	//process remainder of in samples that didn't fill a block
	data.end_of_input = 1;
	data.input_frames  = in.getNumFrames() % blockSize;
	data.src_ratio = 1.0 / safeFactor( in.getTimeOfFrame( (numBlocks + 1) * blockSize ) );
	src_process( state, &data );

	src_delete( state );

	return out;                   
	}
AudioBuffer repitch( ProcInput in, double factor )
	{
	return repitch( in, [factor]( double t ){ return factor; } );
	}

//========================================================
// TODO
//========================================================

//RealFunc ADSR( double A, double D, double S, double R, double Aexp = 0, double Dexp = 0, double Rexp = 0 );

} //end namespace xcdp