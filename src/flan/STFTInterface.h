#pragma once

#include "defines.h"
#include "flan/Utility/MF.h"

namespace flan {

// Interface for time/frequency formats
class STFTInterface
{
public:
	struct Format 
		{ 
		Channel numChannels = 0;
		Frame numFrames = 0;
		Bin numBins = 0;
		SampleRate sampleRate = 48000; 
		Frame hopSize = 0;
		Frame windowSize = 0;
		//window type?
		};

	virtual fFrame timeToFrame( Time ) const = 0;
	virtual Time frameToTime( fFrame ) const = 0;
	virtual fBin frequencyToBin( Frequency ) const = 0;
	virtual Frequency binToFrequency( fBin ) const = 0;

	virtual Channel getNumChannels() const = 0;
	virtual Frame getNumFrames() const = 0;
	virtual Bin	getNumBins( Frame f = -1 ) const = 0;
	virtual SampleRate getSampleRate() const = 0;

	Time getLength() const { return frameToTime( getNumFrames() ); }
	Frequency getHeight() const { return binToFrequency( getNumBins() ); }
	bool isNull() const 
		{ 
		return getSampleRate() > 0
			&& getNumChannels() > 0
			&& getNumFrames() > 0
			&& getNumBins() > 0;
		}
};

}