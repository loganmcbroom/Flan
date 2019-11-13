#pragma once

#include <vector>
#include <string>

#include "Types.h"
#include "AudioBuffer.h"

namespace xcdp {

class PVOC;

class Audio : public AudioBuffer
{
public:
	
	Audio() : AudioBuffer() {} //Only called when erroring out
	Audio( const std::string & filePath ) : AudioBuffer( filePath ) {}
	Audio( const AudioBuffer::Format & other ) : AudioBuffer( other ) {}

	//======================================================
	//	Conversions
	//======================================================

	// Main PVOC analysis function
	PVOC convertToPVOC( const size_t frameSize = 2048, const size_t overlaps = 16 ) const;

	//===========================================================================================
	//	Procs
	//===========================================================================================

	Audio monoToStereo( ) const;

	//Multiply input signal by volumeLevel
	Audio modifyVolume( RealFunc volumeLevel ) const;

	//Normalize input to level
	Audio setVolume( double level = 1.0 ) const;

	//Apply a shaper function to the input samples
	Audio waveshape( RealFunc shaper ) const;

	/** Pan input by panAmount
	 * A panAmount of -1 and 1 correspond to hard left and hard right
	 * This uses sin panning with gain, so a 0 pan won't alter the signal but a non-zero pan can induce clipping
	 */
	Audio pan( RealFunc panAmount ) const;

	/*==============================================================
	 *		Editing
	 *	These processes manipulate discrete pieces of audio in time
	*///============================================================

	//Loops the input n times, applying mod to each loop
	Audio iterate( size_t n, std::function< Audio (const Audio &, size_t n) > mod = 0 ) const;

	// extend freeze
	// drunk walk / scramble
	// reverse
	// extend loop
	// extend repittions? Iterate could call it

	Audio cut( double startTime, double endTime ) const;

	static Audio mix( AudioVec ins, std::vector< RealFunc > balances = std::vector< RealFunc >() );

	//==============================================================

	//==============================================================

	//Scale pitch by factor
	Audio repitch( RealFunc factor ) const;

	Audio convolve( const std::vector<double> & ) const;
};

} // End namespace xcdp
