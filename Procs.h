#pragma once

#include "AudioBuffer.h"

namespace xcdp {

//Multiply input signal by volumeLevel
AudioBuffer modifyVolume( ProcInput in, RealFunc volumeLevel );
AudioBuffer modifyVolume( ProcInput in, double volumeLevel );

//Normalize input to level
AudioBuffer setVolume( ProcInput in, double level = 1.0 );

//Add inputs, each scaled by balances. Balances default to 1/n for n inputs.
AudioBuffer mix( ProcInputVec ins, std::vector< RealFunc > balances = std::vector< RealFunc >() );
AudioBuffer mix( ProcInputVec ins, const std::vector< double > & balances );

//Apply a shaper function to the input samples
AudioBuffer waveshape( ProcInput in, RealFunc shaper );

/** Pan input by panAmount
 * A panAmount of -1 and 1 correspond to hard left and hard right
 * This uses sin panning with gain, so a 0 pan won't alter the signal but a non-zero pan can induce clipping
 */
AudioBuffer pan( ProcInput in, RealFunc panAmount );
AudioBuffer pan( ProcInput in, double panAmount );

//Loops the input n times, applying mod to each loop
AudioBuffer iterate( ProcInput in, int n, std::function< AudioBuffer (ProcInput, int n) > mod = 0 );

//Cut a piece out of the input
AudioBuffer cutAtSamples( ProcInput in, int startSample, int endSample );
AudioBuffer cutAtTimes( ProcInput in, double startTime, double endTime );

//Join any number of input files into one file
AudioBuffer join( AudioBufferVec ins );

//Scale pitch by factor
AudioBuffer repitch( ProcInput in, RealFunc factor );
AudioBuffer repitch( ProcInput in, double factor );

} //end namespace xcdp