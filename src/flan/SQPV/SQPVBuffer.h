#pragma once

#include <utility>
#include <vector>
#include <complex>

#include "flan/defines.h"
#include "flan/Utility/MF.h"

namespace flan {

class SQPVBuffer
{
public:
	struct Format 
		{ 
		Channel numChannels = 0;
		Frame numFrames = 0;
		fBin binsPerOctave = 0;
		SampleRate sampleRate = 48000; 
		std::pair<Frequency, Frequency> bandwidth;
		};

	SQPVBuffer( const SQPVBuffer & ) = delete;
	SQPVBuffer( SQPVBuffer && ) = default;
	SQPVBuffer& operator=( const SQPVBuffer & ) = delete;
	SQPVBuffer& operator=( SQPVBuffer && ) = default;
	~SQPVBuffer() = default;

	SQPVBuffer();
	SQPVBuffer( const Format & );

	MF getMF( Channel, Frame, Bin ) const;
	MF & getMF( Channel, Frame, Bin );

	const Format & getFormat() const;

	fFrame timeToFrame( Time ) const;
	Time frameToTime( fFrame ) const;
	fBin frequencyToBin( Frequency ) const;
	Frequency binToFrequency( fBin ) const;

	Channel getNumChannels() const;
	Frame getNumFrames() const;
	Bin	getNumBins() const;
	SampleRate getSampleRate() const;
	std::pair<Frequency, Frequency> getBandwidth() const;
	fBin getBinsPerOctave() const;
	float getQuality() const;

	Time getLength() const { return frameToTime( getNumFrames() ); }
	bool isNull() const;
	void clearBuffer();
	SQPVBuffer copy() const;

	size_t buffer_access( Channel, Frame, Bin ) const;

private:

	const Format format;
	const Bin numBins;
	const float quality;
	std::vector<MF> buffer;
};

}