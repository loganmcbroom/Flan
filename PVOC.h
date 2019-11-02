#pragma once

#include "Types.h"

namespace xcdp {

class Audio;

class PVOC
{
public:

	//PVOC( const AudioBuffer & );

	//======================================================
	//	I/O
	//======================================================

	//======================================================
	//	Getters
	//======================================================

	size_t getNumChannels() const;
	size_t getNumFrames() const;
	size_t getNumBins() const;
	size_t getFrameSize() const;
	size_t getSampleRate() const;
	double getBinWidth() const;
	double getMaxPartialMagnitude() const;

	//======================================================
	//	Setters
	//======================================================

	void setNumChannels( size_t newNumChannels );
	void setNumFrames( size_t newNumFrames );
	void setNumBins( size_t newNumBins );
	void setBufferSize( size_t newNumChannels,  size_t newNumFrames, size_t newNumBins );
	void setBufferSize( const PVOC & );
	void setFrameSize( size_t newFrameSize );
	void setSampleRate( size_t newSampleRate );
	void copyFormat( const PVOC & );
	void clearBuffer();

	//======================================================
	//	Conversions
	//======================================================

	Audio getAudio() const;
	const PVOC & getSpectrograph( const std::string & fileName = std::string("spectrograph.tga") ) const;

	//======================================================
	//	Procs
	//======================================================

	PVOC timeAverage( int factor ) const;

	PVOC repitch( double factor );
	//PVOC gate( RealFunc cutoff ) const;
	//PVOC gate( double cutoff ) const;
	//PVOC invert( int lowerBound, int upperBound ) const;


	//======================================================
	//	Members
	//======================================================

	size_t sampleRate = 48000;
	size_t overlaps = 4;
	//saving format?
	//window type

	//Main PVOC buffer, stored as:
	//	Channel -> Frame -> Bin -> ( Magnitude, True Frequency )
	struct MFPair { double magnitude, frequency; };
	std::vector< std::vector< std::vector< MFPair > > > buffer;
};

} // End namespace xcdp

