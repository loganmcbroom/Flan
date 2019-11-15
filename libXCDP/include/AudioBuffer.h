#pragma once

#include <vector>
#include <string>

namespace xcdp {

class AudioBuffer
{
public:
	struct Format 
		{
		size_t numChannels = 0, numSamples = 0;
		size_t sampleRate = 48000;
		};
	
	AudioBuffer();
	AudioBuffer( const Format & ); 
	AudioBuffer( const std::string & filename );

	//======================================================
	//	I/O
	//======================================================

	void load( const std::string & filePath );
	bool save( const std::string & filePath ) const;
	void printSummary() const;

	//======================================================
	//	Getters
	//======================================================

	double getSample( size_t channel, size_t frame ) const;
	
	Format getFormat() const;
	size_t getNumChannels() const;
	size_t getNumSamples() const;
	size_t getSampleRate() const;

	double sampleToTime( size_t sample ) const;
	size_t timeToSample( double time ) const;
	double getLength() const;

	double getMaxSampleMagnitude() const;

	//======================================================
	//	Setters
	//======================================================

	void setSample( size_t channel, size_t frame, double sample );
	double & getSample( size_t channel, size_t frame );
	void clearBuffer();

private://=================================================================================================

	size_t getPos( size_t, size_t ) const;

	Format format;
	std::vector< double > buffer;
};

} // End namespace xcdp
