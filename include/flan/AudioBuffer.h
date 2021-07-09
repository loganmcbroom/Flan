#pragma once

#include <vector>
#include <string>
#include <memory>

#include "defines.h"

namespace flan {

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
		Channel numChannels = 0;
		Frame numFrames = 0;
		uint32_t sampleRate = 48000;
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

	/** Returns a deep copy of the AudioBuffer.
	 */
	AudioBuffer deepCopy() const;

	/** Return true if the Audio is in a state which cannot be processed. This includes a 0-channel buffer,
	 *	a 0-frame buffer, or a buffer with a 0 sample rate.
	 */
	bool isNull() const;

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
	Sample getSample( Channel channel, Frame frame ) const;

	/** Buffer-safe sample access method.
	 *	\param channel
	 *	\param frame
	 */
	Sample getSample_s( Channel channel, Frame frame ) const;
	
	/** Returns the format being used.
	 */
	Format getFormat() const;

	/** Returns the number of channels.
	 */
	Channel getNumChannels() const;

	/** Returns the number of frames.
	 */
	Frame getNumFrames() const;

	/** Returns the sample rate.
	 */
	uint32_t getSampleRate() const;

	/** Returns length in seconds.
	 */
	Time getLength() const;

	/** Returns the magnitude of whichever sample has the largest magnitude in the given range.
	 *	\param startTime Start of range.
	 *	\param endTime End of range. Passing 0 will use the maximum time.
	 */
	Sample getMaxSampleMagnitude( Time startTime = 0, Time endTime = 0 ) const;

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
	void setSample( Channel channel, Frame frame, Sample sample );

	/** Additional sample setter method. It is occasionally syntactically easier to assign samples by reference access.
	 *	\param channel
	 *	\param frame
	 */
	Sample & getSample( Channel channel, Frame frame );
	
	/** This zero fills the buffer.
	 */
	void clearBuffer();

	/** Raw buffer access. Inderect access is preferred, but for computationally expensive 
	 *	transformations, raw access usually provides an optimization.
	 *	\param channel
	 *	\param frame
	 */
	Sample * getSamplePointer( Channel channel, Frame frame );

	/** Read-only raw buffer access. Inderect access is preferred, but for computationally expensive 
	 *	transformations, raw access usually provides an optimization.
	 *	\param channel
	 *	\param frame
	 */
	const Sample * getSamplePointer( Channel channel, Frame frame ) const;

	std::vector<Sample>::iterator begin();
	std::vector<Sample>::iterator end();
	std::vector<Sample>::const_iterator begin() const;
	std::vector<Sample>::const_iterator end() const;

private://=================================================================================================

	inline uint32_t getBufferPos( Channel, Frame ) const;

	Format format;
	std::shared_ptr<std::vector<Sample>> buffer;
};

/** Serialization.
 */
std::ostream & operator<<( std::ostream & os, const AudioBuffer & audio );

} // End namespace flan
