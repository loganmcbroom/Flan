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
	
	Audio( const std::string & filePath ) : AudioBuffer( filePath ) {}
	Audio( const AudioBuffer::Format & other ) : AudioBuffer( other ) {}

	//======================================================
	//	Conversions
	//======================================================

	PVOC convertToPVOC( const size_t frameSize = 2048, const size_t overlaps = 16 ) const;

	//======================================================
	//	Procs
	//======================================================

	Audio mono2Stereo( ) const;

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

	//Loops the input n times, applying mod to each loop
	Audio iterate( size_t n, std::function< Audio (const Audio &, size_t n) > mod = 0 ) const;

	//Cut a piece out of the input
	Audio cutAtSamples( size_t startSample, size_t endSample ) const;
	Audio cutAtTimes( double startTime, double endTime ) const;

	//Scale pitch by factor
	Audio repitch( RealFunc factor ) const;

	Audio convolve( const std::vector<double> & ) const;
};

} // End namespace xcdp
