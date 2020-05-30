#pragma once

#include <vector>
#include <string>
#include <memory>

namespace xcdp {

/** AudioBuffer stores audio data and provides basic buffer access, conversion constants, loading, and saving.
 *	Access to the raw sample buffer is given, but AudioBuffer::getSample and 
 *	AudioBuffer::setSample are preferred when speed is not a factor.
 *	Data is stored in channel-major order, in other words, the entire first channel 
 *	is stored in memory before the second channel's storage begins.
 *
 *	A note on terminology, the term "frame" is used to describe points in time at which samples are taken.
 *	The term "sample" is used to describe the actual data at that point in time.
 */
class AudioBuffer
{
public:
	/** AudioBuffer::Format stores all required information outside of the buffer. 
	 *	This is used to easily transfer the way an AudioBuffer is stored on disk to a new AudioBuffer without copying the buffer.
	 *	In practice, a function transforming an AudioBuffer will call AudioBuffer::getFormat(), modify members as needed, and
	 *	then construct a blank AudioBuffer with the Format constructor.
	 */
	struct Format 
		{
		size_t numChannels = 0;
		size_t numFrames = 0;
		size_t sampleRate = 48000;
		};
	
	/** Default constructor
	 */
	AudioBuffer();

	/** Format copy constructor, constructs an AudioBuffer with the given format and an uninitialized buffer.
	 *	\param format Format to use in construction.
	 */
	AudioBuffer( const Format & format ); 

	/** Load constructor, calls AudioBuffer::load.
	 *	\param filePath File to load. Accepts any format accepted by libsndfile. Notably cannot load mp3.
	 */
	AudioBuffer( const std::string & filePath );

	//======================================================
	//	I/O
	//======================================================

	/** File loading.
	 *	\param filePath File to load. If libsndfile is enabled this can be any format accepted by libsndfile, notably not mp3.
	 *		If libsndfile is disabled this can only be a wave file.
	 */
	bool load( const std::string & filePath );

	/** File saving.
	 *	\param filePath File path to save at.
	 *	\param format The libsndfile format to save to. The default of -1 will save as 24bit PCM WAVE.
	 *		If libsndfile is disabled this does nothing, and 24bit PCM WAVE will always be used.
	 */
	bool save( const std::string & filePath, int format = -1 ) const;

	/** Prints buffer dimensions and sample rate to cout.
	 */
	void printSummary() const;

	//======================================================
	//	Getters
	//======================================================

	/** Main sample access method.
	 *	\param channel
	 *	\param frame
	 */
	float getSample( size_t channel, size_t frame ) const;
	
	/** Returns the format being used.
	 */
	Format getFormat() const;

	/** Returns the number of channels.
	 */
	size_t getNumChannels() const;

	/** Returns the number of frames.
	 */
	size_t getNumFrames() const;

	/** Returns the sample rate.
	 */
	size_t getSampleRate() const;

	/** Returns length in seconds.
	 */
	float getLength() const;

	/** Returns the magnitude of whichever sample has the largest magnitude.
	 */
	float getMaxSampleMagnitude() const;

	/** Returns a unit fraction for converting frames to seconds.
	 */
	float frameToTime() const;

	/** Returns a unit fraction for converting seconds to frames
	 */
	float timeToFrame() const;


	//======================================================
	//	Setters
	//======================================================

	/** Main sample setter method.
	 *	\param channel
	 *	\param frame
	 *	\param sample The value to write.
	 */
	void setSample( size_t channel, size_t frame, float sample );

	/** Additional sample setter method. It is occasionally syntactically easier to assign samples by reference access.
	 *	\param channel
	 *	\param frame
	 */
	float & getSample( size_t channel, size_t frame );
	
	/** This zero fills the buffer.
	 */
	void clearBuffer();

	/** Raw buffer access. Inderect access is preferred, but for computationally expensive 
	 *	transformations, raw access usually provides an optimization.
	 *	\param channel
	 *	\param frame
	 */
	float * getSamplePointer( size_t channel, size_t frame );

	/** Read-only raw buffer access. Inderect access is preferred, but for computationally expensive 
	 *	transformations, raw access usually provides an optimization.
	 *	\param channel
	 *	\param frame
	 */
	const float * getSamplePointer( size_t channel, size_t frame ) const;

private://=================================================================================================

	inline size_t getBufferPos( size_t, size_t ) const;

	Format format;
	std::shared_ptr<std::vector<float>> buffer;
};

/** Serialization.
 */
std::ostream & operator<<( std::ostream & os, const AudioBuffer & audio );

} // End namespace xcdp
