#pragma once

#include <vector>
#include <string>
#include <memory>

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

	float getSample( size_t channel, size_t frame ) const;
	
	Format getFormat() const;
	size_t getNumChannels() const;
	size_t getNumSamples() const;
	size_t getSampleRate() const;

	float sampleToTime( size_t sample ) const;
	size_t timeToSample( float time ) const;
	float getLength() const;

	float getMaxSampleMagnitude() const;

	//======================================================
	//	Setters
	//======================================================

	void setSample( size_t channel, size_t frame, float sample );
	float & getSample( size_t channel, size_t frame );
	
	void clearBuffer();

	//Raw buffer access methods, only use if it adds significant speed
	std::shared_ptr<std::vector<float>> getBuffer() { return buffer; }
	const std::shared_ptr<std::vector<float>> getBuffer() const { return buffer; }
	float * getSamplePointer( size_t channel, size_t frame );
	const float * getSamplePointer( size_t channel, size_t frame ) const;

private://=================================================================================================

	inline size_t getPos( size_t, size_t ) const;

	Format format;
	std::shared_ptr<std::vector<float>> buffer;
};

} // End namespace xcdp
