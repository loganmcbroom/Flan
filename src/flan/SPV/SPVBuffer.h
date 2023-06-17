#pragma once

#include <vector>

#include "flan/Utility/MF.h"

namespace flan {

class Audio;

class SPVBuffer
{
public:
	struct Format 
		{ 
		Channel numChannels = 0;
		Frame numFrames = 0;
		Bin numBins = 0;
		FrameRate sampleRate = 48000;
		};

	SPVBuffer( const SPVBuffer & ) = delete;
	SPVBuffer( SPVBuffer && ) = default;
	SPVBuffer& operator=( const SPVBuffer & ) = delete;
	SPVBuffer& operator=( SPVBuffer && ) = default;
	~SPVBuffer() = default;

	SPVBuffer();
	SPVBuffer( Format );

	const Format & getFormat() const;
	fFrame timeToFrame( Second ) const;
	Second frameToTime( fFrame ) const;
	fBin frequencyToBin( Frequency ) const;
	Frequency binToFrequency( fBin ) const;

	Channel getNumChannels() const;
	Frame getNumFrames() const;
	Bin	getNumBins() const;
	FrameRate getSampleRate() const;
	FrameRate getAnalysisRate() const;

	Second getLength() const { return frameToTime( getNumFrames() ); }
	Frequency getHeight() const { return binToFrequency( getNumBins() ); }
	bool isNull() const;

	MF getMF( Channel channel, Frame frame, Bin bin ) const;
	MF & getMF( Channel channel, Frame frame, Bin bin );
	void clearBuffer();
	SPVBuffer copy() const;

	size_t getBufferPos( Channel, Frame, Bin ) const;

	// Avoid using this when possible
	std::vector<MF> & getBuffer();

private:

	Format format;
	std::vector<MF> buffer;
};

}