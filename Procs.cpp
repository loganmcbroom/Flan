#include "Procs.h"

#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>

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

//Copy sample rate and bit depth of input
AudioFile<double> copyAudioSpec( ProcInput in )
	{
	AudioFile<double> out;
	out.setSampleRate( in.getSampleRate() );
	out.setBitDepth( in.getBitDepth() );
	return out;
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
		return a.getNumSamplesPerChannel() < b.getNumSamplesPerChannel();
		} )->getNumSamplesPerChannel();
	}

double getMaxSample( ProcInput in )
	{
	double maxSample = 0.0;
	
	for( int channel = 0; channel < in.getNumChannels(); ++channel )
		{
		double maxSampleInChannel = std::abs( *std::max_element( 
			in.samples[channel].begin(), 
			in.samples[channel].end(), 
			[](double a, double b){ return std::abs(a) < std::abs(b); }));
		maxSample = std::max( maxSample, maxSampleInChannel );
		}

	return maxSample;
	}

//========================================================
// Procs
//========================================================

AudioFile<double> modifyVolume( ProcInput in, RealFunc volumeLevel )
	{
	auto out = copyAudioSpec( in );
	out.setAudioBufferSize( in.getNumChannels(), in.getNumSamplesPerChannel() );

	for( int channel = 0; channel < in.getNumChannels(); ++channel )
		for( int sample = 0; sample < in.getNumSamplesPerChannel(); ++sample )
			{
			double calculatedSample = in.samples[channel][sample]*volumeLevel(sample);
			out.samples[channel][sample] = std::clamp( calculatedSample, -1.0, 1.0 );
			}

	return out;
	}
AudioFile<double> modifyVolume( ProcInput in, double volumeLevel )
	{
	return modifyVolume( in, [ volumeLevel ]( double ){ return volumeLevel; } );
	}

AudioFile<double> setVolume( ProcInput in, double level )
	{
	// Divide by maxVolume to normalize, multiply by level to set
	double volume = level / getMaxSample( in );
	return modifyVolume( in, [ volume ]( double t ){ return volume; } );
	}

AudioFile<double> mix( ProcInputVec ins, std::vector< RealFunc > balances )
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
		printToLog( "FIX: You can't mix files with differing channel counts at the moment. \n" );
		return AudioFile<double>();
		}

	int numChannels = getMaxNumChannels( ins );
	int numSamples  = getMaxNumSamples ( ins );

	auto out = copyAudioSpec( ins[0] );
	out.setAudioBufferSize( numChannels, numSamples );

	for( int channel = 0; channel < numChannels; ++channel )
		for( int sample = 0; sample < numSamples; ++sample )
			{
			double mixedSample = 0.0;
			for( int i = 0; i < ins.size(); ++i )
				mixedSample += balances[i]( sample ) * ins[i].samples[channel][sample];
			out.samples[channel][sample] = std::clamp( mixedSample, -1.0, 1.0 );
			}

	return out;
	}
AudioFile<double> mix( ProcInputVec ins, const std::vector< double > & balances )
	{
	std::vector< RealFunc > funcs( balances.size() );
	for( unsigned int i = 0; i < balances.size(); ++i )
		funcs[i] = [i, &balances](double){ return balances[i]; };

	return mix( ins, funcs );
	}

AudioFile<double> waveshape( ProcInput in, RealFunc shaper )
	{
	auto out = copyAudioSpec( in );
	out.setAudioBufferSize( in.getNumChannels(), in.getNumSamplesPerChannel() );

	for( int channel = 0; channel < in.getNumChannels(); ++channel )
		for( int sample = 0; sample < in.getNumSamplesPerChannel(); ++sample )
			{
			double calculatedSample = shaper( in.samples[channel][sample] );
			out.samples[channel][sample] = std::clamp( calculatedSample, -1.0, 1.0 );
			}

	return out;
	}

//Duplicates a mono channel as two stereo channels
AudioFile<double> mono2Stereo( ProcInput in )
	{
	if( in.getNumChannels() != 1 )
		printToLog( "That wasn't a mono file you tried to make stereo just now.\n" );

	auto out = in;
	out.setNumChannels( 2 );
	out.samples[2] = out.samples[1];
	return out;
	}

AudioFile<double> pan( ProcInput in, RealFunc panAmount )
	{
	//Stereo panning algorithm
	auto stereoPan = []( ProcInput in, RealFunc panAmount )
		{
		auto out = copyAudioSpec( in );
		out.setAudioBufferSize( in.getNumChannels(), in.getNumSamplesPerChannel() );

		for( int channel = 0; channel < in.getNumChannels(); ++channel )
			for( int sample = 0; sample < in.getNumSamplesPerChannel(); ++sample )
				{
				double calculatedSample = in.samples[channel][sample]
					* sin( pi / 4.0 * ( std::clamp( panAmount(sample), -1.0, 1.0 ) + 3.0 - double(channel)*2.0 ) )
					* sqrt( 2.0 );
				out.samples[channel][sample] = std::clamp( calculatedSample, -1.0, 1.0 );
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
AudioFile<double> pan( ProcInput in, double panAmount )
	{
	return pan( in, [panAmount](double){ return panAmount; } );
	}

AudioFile<double> loop( ProcInput in, int n )
	{
	auto out = copyAudioSpec( in );
	out.setNumChannels( in.getNumChannels() );
	out.setNumSamplesPerChannel( in.getNumSamplesPerChannel() * n );

	for( int channel = 0; channel < out.getNumChannels(); ++channel )
		for( int sample = 0; sample < out.getNumSamplesPerChannel(); ++sample )
			out.samples[channel][sample] = in.samples[channel][sample % in.getNumSamplesPerChannel()];

	return out;
	}

AudioFile<double> cutAtSamples( ProcInput in, int startSample, int endSample )
	{
	//input validity checking
	startSample = std::min( 0, startSample );
	endSample = std::max( in.getNumSamplesPerChannel(), endSample );
	if( endSample < startSample ) endSample = startSample;

	auto out = copyAudioSpec( in );
	out.setAudioBufferSize( in.getNumChannels(), endSample - startSample );

	for( int channel = 0; channel < out.getNumChannels(); ++channel )
		for( int sample = 0; sample < out.getNumSamplesPerChannel(); ++sample )
			out.samples[channel][sample] = in.samples[channel][ startSample + sample ];

	return out;
	}
AudioFile<double> cutAtTimes( ProcInput in, double startTime, double endTime )
	{
	//Input time validity is checked in cutAtSamples
	return cutAtSamples( in, startTime / in.getSampleRate(), endTime / in.getSampleRate() );
	}

AudioFile<double> join( AudioFileVec ins )
	{
	int numChannels = getMaxNumChannels( ins );
	int numSamples = 0;
	for( auto & in : ins ) numSamples += in.getNumSamplesPerChannel();

	auto out = copyAudioSpec( ins[0] );
	out.setAudioBufferSize( numChannels, numSamples );

	int currentOutStartSample = 0;

	//For each in, copy in into out
	for( int in = 0; in < ins.size(); ++in )
		{
		for( int channel = 0; channel < ins[in].getNumChannels(); ++channel )
			for( int sample = 0; sample < ins[in].getNumSamplesPerChannel(); ++sample )
				out.samples[channel][currentOutStartSample + sample] = ins[in].samples[channel][sample];
		currentOutStartSample += ins[in].getNumSamplesPerChannel();
		}

	return out;
	}