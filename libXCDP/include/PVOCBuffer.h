#pragma once

#include <vector>
#include <string>
#include <memory>

namespace xcdp {

class PVOCBuffer
{
public:
	struct MFPair { float magnitude, frequency; };
	struct Format 
		{ 
		size_t numChannels = 0, numFrames = 0, numBins = 0;
		size_t sampleRate = 48000; 
		size_t overlaps = 32;
		//saving format?
		//window type?
		};

	PVOCBuffer( const Format & );
	PVOCBuffer( const std::string & filename );

	//======================================================
	//	I/O
	//======================================================

	bool load( const std::string & filename );
	bool save( const std::string & filename ) const;

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
	float getBinWidth() const;
	float getMaxPartialMagnitude() const;

	size_t timeToFrame( float time ) const;
	float timeToFrame_d( float time ) const;
	float frameToTime( size_t frame ) const;

	size_t frequencyToBin( float freq ) const;
	float frequencyToBin_d( float freq ) const;
	float binToFrequency( size_t bin ) const;

	//======================================================
	//	Setters
	//======================================================

	void setBin( size_t channel, size_t frame, size_t bin, MFPair value );
	MFPair & getBin( size_t channel, size_t frame, size_t bin );
	void clearBuffer();

	//raw buffer access, use with caution
	std::shared_ptr<std::vector< MFPair >> getBuffer() { return buffer; }
	const std::shared_ptr<std::vector< MFPair >> getBuffer() const { return buffer; }

private://=================================================================================================

	size_t getPos( size_t, size_t, size_t ) const;
	
	Format format;
	std::shared_ptr<std::vector< MFPair >> buffer;
};

} // End namespace xcdp

