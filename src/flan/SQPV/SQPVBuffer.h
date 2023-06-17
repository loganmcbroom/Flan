#pragma once

#include <utility>
#include <vector>
#include <complex>

#include "flan/defines.h"
#include "flan/Utility/MP.h"

namespace flan {

class SQPVBuffer
{
public:
	struct Format 
		{ 
		Channel numChannels = 0;
		Frame numFrames = 0;
		fBin binsPerOctave = 0;
		FrameRate sampleRate = 48000; 
		std::pair<Frequency, Frequency> bandwidth;
		};

	SQPVBuffer( const SQPVBuffer & ) = delete;
	SQPVBuffer( SQPVBuffer && ) = default;
	SQPVBuffer& operator=( const SQPVBuffer & ) = delete;
	SQPVBuffer& operator=( SQPVBuffer && ) = default;
	~SQPVBuffer() = default;

	SQPVBuffer();
	SQPVBuffer( const Format & );

	MP getMP( Channel, Frame, Bin ) const;
	MP & getMP( Channel, Frame, Bin );

	const Format & getFormat() const;

	fFrame timeToFrame( Second ) const;
	Second frameToTime( fFrame ) const;
	fBin frequencyToBin( Frequency ) const;
	Frequency binToFrequency( fBin ) const;
	Pitch frequencyToPitch( Frequency ) const;
	Frequency pitchToFrequency( Pitch ) const;
	UnsignedPitch binToPitch( fBin ) const;
	fBin pitchToBin( UnsignedPitch ) const;

	Channel getNumChannels() const;
	Frame getNumFrames() const;
	Bin	getNumBins() const;
	FrameRate getSampleRate() const;
	FrameRate getAnalysisRate() const;
	std::pair<Frequency, Frequency> getFrequencyBandwidth() const;
	std::pair<UnsignedPitch, UnsignedPitch> getPitchBandwidth() const;
	fBin getBinsPerOctave() const;
	Cycle getQ() const;
	Frequency getBinFrequency( Bin ) const;
	Second getLength() const { return frameToTime( getNumFrames() ); }
	bool isNull() const;
	void clearBuffer();
	SQPVBuffer copy() const;
	Magnitude getMaxPartialMagnitude() const;
	Magnitude getMaxPartialMagnitude( Frame startFrame, Frame endFrame = 0, Bin startBin = 0, Bin endBin = 0 ) const;
	Frame getPeriod( Bin ) const;

	size_t getBufferPos( Channel, Frame, Bin ) const;

private:

	const Format format;

	// Commonly used data is computed once at construction
	const std::pair<UnsignedPitch, UnsignedPitch> pitch_bandwidth;
	const Bin numBins;
	const Cycle Q; // "The number of cycles needed to make an analysis" - Sliding With A Constant Q
	const std::vector<Frequency> bin_frequencies;

	std::vector<MP> buffer;
};

}