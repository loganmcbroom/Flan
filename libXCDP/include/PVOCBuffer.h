#pragma once

#include <vector>

namespace xcdp {

class PVOCBuffer
{
public:
	struct MFPair { double magnitude, frequency; };
	struct Format 
		{ 
		size_t numChannels = 0, numFrames = 0, numBins = 0;
		size_t sampleRate = 48000; 
		size_t overlaps = 32;
		//saving format?
		//window type?
		};

	PVOCBuffer( const Format & );

	//======================================================
	//	I/O
	//======================================================

	//======================================================
	//	Getters
	//======================================================

	MFPair getBin( size_t channel, size_t frame, size_t bin ) const;

	Format getFormat() const;
	size_t getNumChannels() const;
	size_t getNumFrames() const;
	size_t getNumBins() const;
	size_t getSampleRate() const;
	size_t getOverlaps() const;

	size_t getFrameSize() const;
	double getBinWidth() const;
	double getMaxPartialMagnitude() const;
	long long timeToFrame( double t ) const;
	double frameToTime( size_t frame ) const;

	//======================================================
	//	Setters
	//======================================================

	void setBin( size_t channel, size_t frame, size_t bin, MFPair value );
	void clearBuffer();

private://=================================================================================================

	size_t getPos( size_t, size_t, size_t ) const;
	
	const Format format;
	std::vector< MFPair > buffer;
};

} // End namespace xcdp

