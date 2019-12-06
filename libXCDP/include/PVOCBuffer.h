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
	PVOCBuffer getFrame( size_t frame ) const;

	Format getFormat() const;
	size_t getNumChannels() const;
	size_t getNumFrames() const;
	size_t getNumBins() const;
	size_t getSampleRate() const;
	size_t getOverlaps() const;

	size_t getFrameSize() const;
	double getBinWidth() const;
	double getMaxPartialMagnitude() const;

	size_t timeToFrame( double time ) const;
	double timeToFrame_d( double time ) const;
	double frameToTime( size_t frame ) const;

	size_t frequencyToBin( double freq ) const;
	double frequencyToBin_d( double freq ) const;
	double binToFrequency( size_t bin ) const;

	//======================================================
	//	Setters
	//======================================================

	void setBin( size_t channel, size_t frame, size_t bin, MFPair value );
	MFPair & getBin( size_t channel, size_t frame, size_t bin );
	void clearBuffer();

private://=================================================================================================

	size_t getPos( size_t, size_t, size_t ) const;
	
	const Format format;
	std::vector< MFPair > buffer;
};

} // End namespace xcdp

