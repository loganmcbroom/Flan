#pragma once

#include <vector>
#include <string>

#include <sndfile.h>

#include "Types.h"

namespace xcdp 
{

class SpectralBuffer;

class AudioBuffer
{
public:
	AudioBuffer();
	AudioBuffer( const std::string & filePath );
	AudioBuffer( const SpectralBuffer & spectra );

	struct ChannelProxy;
	const ChannelProxy operator[]( int channel ) const;


	//======================================================
	//	I/O
	//======================================================

	//Load a file into the buffer
	bool load( const std::string & filePath );

	//Save the buffer into a file. Non-const to allow clipping on save.
	bool save( const std::string & filePath );

	//Print some buffer data to the console
	void printSummary() const;



	//======================================================
	//	Getters
	//======================================================

	//Get the sample value at the specified frame and channel
	double getSample( int channel, int frame ) const;

	//Get the current number of channels
	int getNumChannels() const;

	//Get the current number of frames
	int getNumFrames() const;

	//Get the sample rate
	int getSampleRate() const;

	//Return the time at which the sample occurs
	double getTimeOfFrame( int frame ) const;

	//Returns a single channel of the buffer (WARNING: This is expensive)
	std::vector<double> getChannel( int channel );

	//Get the spectrum via dft
	SpectralBuffer getSpectral() const;



	//======================================================
	//	Setters
	//======================================================

	//Set the sample value at the specified frame and channel
	void setSample( int channel, int frame, double sample );

	//Set the current number of channels. 
	// This can take a long time if data is already loaded in.
	void setNumChannels( int numChannels );

	//Set the current number of frames
	void setNumFrames( int numFrames );

	//Set the current channel count and frame count
	void setBufferSize( int numChannels, int numFrames );

	//Set the sample rate
	void setSampleRate( int sampleRate );

	//Copy the sample rate and file format of another AudioBuffer
	void copyFormat( const AudioBuffer & other );



	//======================================================
	//	Members
	//======================================================

	//Metadata struct
	SF_INFO info;

	//Audio data stored in frame major order, i.e. this data:
	// Channel 1 | 11 12 13 14
	// Channel 2 | 21 22 23 24
	// Is stored like this:
	// 11 21 12 22 13 23 14 24
	//If you modify this, as you sometimes need to do, you must also update the info struct
	std::vector<double> buffer;

	//Syntactic sugar array access
	struct ChannelProxy 
		{
		ChannelProxy( const AudioBuffer & _buddy, int _channel ) : buddy( _buddy ), channel( _channel ) {}
		double operator[]( int bin ) const { return buddy.getSample( channel, bin ); };

		const AudioBuffer & buddy;
		const int channel;
		};
};

}