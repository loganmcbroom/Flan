#pragma once

#include <functional>

#include "AudioFile.h"

typedef std::function< double ( int sample ) > RealFunc;
typedef std::vector<AudioFile<double>> AudioFileVec;
typedef const AudioFile<double> & ProcInput;
typedef const AudioFileVec & ProcInputVec;

//Returns the size of the largest sample
double getMaxSample( ProcInput in );

//Multiply input signal by volumeLevel
AudioFile<double> modifyVolume( ProcInput in, RealFunc volumeLevel );
AudioFile<double> modifyVolume( ProcInput in, double volumeLevel );

//Normalize input to level
AudioFile<double> setVolume( ProcInput in, double level = 1.0 );

//Add inputs, each scaled by balances. Balances default to 1/n for n inputs.
AudioFile<double> mix( ProcInputVec ins, std::vector< RealFunc > balances = std::vector< RealFunc >() );
AudioFile<double> mix( ProcInputVec ins, const std::vector< double > & balances );

//Apply a shaper function to the input samples
AudioFile<double> waveshape( ProcInput in, RealFunc shaper );

/** Pan input by panAmount
 * A panAmount of -1 and 1 correspond to hard left and hard right
 * This uses sin panning with gain, so a 0 pan won't alter the signal but a non-zero pan can induce clipping
 */
AudioFile<double> pan( ProcInput in, RealFunc panAmount );
AudioFile<double> pan( ProcInput in, double panAmount );

//Loops the input n times
AudioFile<double> loop( ProcInput in, int n );

//Cut a piece out of the input
AudioFile<double> cutAtSamples( ProcInput in, int startSample, int endSample );
AudioFile<double> cutAtTimes( ProcInput in, double startTime, double endTime );

//Join any number of input files into one file
AudioFile<double> join( AudioFileVec ins );