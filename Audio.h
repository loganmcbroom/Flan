#pragma once

#include <vector>
#include <string>

#include "Types.h"

/* Figure these bad boys out

//Join any number of input files into one file
AudioBuffer join( AudioBufferVec ins );

//Add inputs, each scaled by balances. Balances default to 1/n for n inputs.
AudioBuffer mix( ProcInputVec ins, std::vector< RealFunc > balances = std::vector< RealFunc >() );
AudioBuffer mix( ProcInputVec ins, const std::vector< double > & balances );

*/

namespace xcdp {

class PVOC;

class Audio
{
public:
	
	Audio();
	Audio( const std::string & filePath );
	//Audio( const PVOC & spectra );

	//Syntactic sugar array access, needs work
	struct ChannelProxy 
		{
		ChannelProxy( const Audio & _buddy, size_t _channel ) : buddy( _buddy ), channel( _channel ) {}
		double operator[]( size_t bin ) const { return buddy.getSample( channel, bin ); };

		const Audio & buddy;
		const size_t channel;
		};
	//const ChannelProxy operator[]( size_t channel ) const;


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
	double getSample( size_t channel, size_t frame ) const;

	//Get the current number of channels
	size_t getNumChannels() const;

	//Get the current number of frames
	size_t getNumFrames() const;

	//Get the sample rate
	size_t getSampleRate() const;

	//Return the time at which the sample occurs
	double getTimeOfFrame( size_t frame ) const;

	//Returns a single channel of the buffer (WARNING: This is expensive)
	std::vector<double> getChannel( size_t channel );

	const std::vector<std::vector<double>> & getBuffer() const;


	//======================================================
	//	Setters
	//======================================================

	//Set the sample value at the specified frame and channel
	void setSample( size_t channel, size_t frame, double sample );

	//Set the current number of channels. 
	void setNumChannels( size_t numChannels );

	//Set the current number of frames
	void setNumFrames( size_t numFrames );

	//Set the current channel count and frame count
	void setBufferSize( size_t numChannels, size_t numFrames );

	//Set the current channel count and frame count
	void setBufferSize( const Audio & other );

	//Set the sample rate
	void setSampleRate( size_t sampleRate );

	//Copy the sample rate and file format of another Audio
	void copyFormat( const Audio & other );

	//======================================================
	//	Conversions
	//======================================================

	PVOC getPVOC( const size_t frameSize = 1024, const size_t overlaps = 4 ) const;

	//======================================================
	//	Procs
	//======================================================

	double getMaxSample() const;
	Audio mono2Stereo( ) const;

	//Multiply input signal by volumeLevel
	Audio modifyVolume( RealFunc volumeLevel ) const;
	Audio modifyVolume( double volumeLevel ) const;

	//Normalize input to level
	Audio setVolume( double level = 1.0 ) const;

	//Apply a shaper function to the input samples
	Audio waveshape( RealFunc shaper ) const;

	/** Pan input by panAmount
	 * A panAmount of -1 and 1 correspond to hard left and hard right
	 * This uses sin panning with gain, so a 0 pan won't alter the signal but a non-zero pan can induce clipping
	 */
	Audio pan( RealFunc panAmount ) const;
	Audio pan( double panAmount ) const;

	//Loops the input n times, applying mod to each loop
	Audio iterate( size_t n, std::function< Audio (const Audio &, size_t n) > mod = 0 ) const;

	//Cut a piece out of the input
	Audio cutAtSamples( size_t startSample, size_t endSample ) const;
	Audio cutAtTimes( double startTime, double endTime ) const;

	//Scale pitch by factor
	Audio repitch( RealFunc factor ) const;
	Audio repitch( double factor ) const;

	//======================================================
	//	Members
	//======================================================

	size_t sampleRate;
	//format?

	//Audio data buffer stored as: Channel -> Frame
	std::vector<std::vector<double>> buffer;

};

} // End namespace xcdp
