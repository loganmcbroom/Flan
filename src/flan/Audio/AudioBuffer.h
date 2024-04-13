#pragma once

#include <vector>
#include <string>
#include <memory>

#include "flan/defines.h"

namespace flan {

/** AudioBuffer stores audio data and provides basic buffer access, conversion constants, loading, and saving.
 *	Access to the raw sample buffer is given, but AudioBuffer::get_sample and 
 *	AudioBuffer::set_sample are preferred when speed is not a factor.
 *	Data is stored in channel-major order, in other words, the entire first channel 
 *	is stored in memory before the second channel's storage begins.
 *
 *	A note on terminology, the term "frame" is used to describe points in time at which samples are taken.
 *	The term "sample" is used to describe the actual data at that point in time.
 */
class AudioBuffer
{
public:
	AudioBuffer( const AudioBuffer & ) = delete;
	AudioBuffer( AudioBuffer && ) = default;
	AudioBuffer& operator=( const AudioBuffer & ) = delete;
	AudioBuffer& operator=( AudioBuffer && ) = default;
	~AudioBuffer() = default;

	/** AudioBuffer::Format stores all required information outside of the buffer. 
	 *	This is used to easily transfer the way an AudioBuffer is stored on disk to a new AudioBuffer without copying the buffer.
	 *	In practice, a function transforming an AudioBuffer will call AudioBuffer::get_format(), modify members as needed, and
	 *	then construct a blank AudioBuffer with the Format constructor.
	 */
	struct Format 
		{
		Channel num_channels = 0;
		Frame num_frames = 0;
		FrameRate sample_rate = 48000;
		};

	struct SndfileStrings 
		{
		std::string title = "";
		std::string copyright = "";
		std::string software = "";
		std::string artist = "";
		std::string comment = "";
		std::string date = "";
		std::string album = "";
		std::string license = "";
		std::string tracknumber = "";
		std::string genre = "";
		};

	struct ChannelSlice
		{
		using Iter = std::vector<Sample>::iterator;

		Iter begin() { return m_begin; }
		Iter end() { return m_end; }

		Iter m_begin;
		Iter m_end;
		};
	
	/** Default constructor
	 */
	AudioBuffer();

	/** Construct from temporary buffer.
	 */
	AudioBuffer( std::vector<float> && buffer, Channel num_channels, FrameRate = 48000 );

	/** Format copy constructor, constructs an AudioBuffer with the given format and an uninitialized buffer.
	 *	\param format Format to use in construction.
	 */
	AudioBuffer( const Format & format ); 

	/** Load constructor, calls AudioBuffer::load.
	 *	\param filepath File to load. Accepts any format accepted by libsndfile. Notably cannot load mp3.
	 */
	AudioBuffer( const std::string & filepath );
	AudioBuffer( const std::string & filepath, SndfileStrings & );

	/** Returns a deep copy of the AudioBuffer.
	 */
	AudioBuffer copy() const;

	/** Return true if the Audio is in a state which cannot be processed. This includes a 0-channel buffer,
	 *	a 0-frame buffer, or a buffer with a 0 sample rate.
	 */
	bool is_null() const;

	/// @brief Scans all samples for nan or inf.
	/// @return Returns true if any sample is nan or inf.
	bool is_nan_or_inf() const;

	//======================================================
	//	I/O
	//======================================================

	/** File loading.
	 *	\param filepath File to load. If libsndfile is enabled this can be any format accepted by libsndfile, notably not mp3.
	 *		If libsndfile is disabled this can only be a wave file.
	 */
	bool load( 
		const std::string & filepath
		);

	bool load( 
		const std::string & filepath, 
		SndfileStrings &
		);

	/** File saving.
	 *	\param filepath File path to save at.
	 *	\param format The libsndfile format to save to. The default of -1 will save as 24bit PCM WAVE.
	 *		If libsndfile is disabled this does nothing, and 24bit PCM WAVE will always be used.
	 */
	bool save( 
		const std::string & filepath, 
		int format = -1, 
		SndfileStrings = SndfileStrings() // Used for string smuggling
		) const;

	/** Prints buffer dimensions and sample rate to cout.
	 */
	void print_summary() const;

	//======================================================
	//	Getters
	//======================================================

	/** Main sample access method.
	 *	\param channel
	 *	\param frame
	 */
	Sample get_sample( Channel channel, Frame frame ) const;
	
	/** Returns the format being used.
	 */
	Format get_format() const;

	/** Returns the number of channels.
	 */
	Channel get_num_channels() const;

	/** Returns the number of frames.
	 */
	Frame get_num_frames() const;

	/** Returns the sample rate.
	 */
	FrameRate get_sample_rate() const;

	/** Returns length in seconds.
	 */
	Second get_length() const;

	/** Returns the magnitude of whichever sample has the largest magnitude in the given range.
	 *	\param start_time Start of range.
	 *	\param end_time End of range. Passing 0 will use the maximum time.
	 */
	Sample get_max_sample_magnitude( Second start_time = 0, Second end_time = 0 ) const;

	/** Returns a unit fraction for converting frames to seconds.
	 */
	Second frame_to_time( fFrame ) const;

	/** Returns a unit fraction for converting seconds to frames
	 */
	fFrame time_to_frame( Second ) const;


	//======================================================
	//	Setters
	//======================================================

	/** Main sample setter method.
	 *	\param channel
	 *	\param frame
	 *	\param sample The value to write.
	 */
	void set_sample( Channel channel, Frame frame, Sample sample );

	/** Additional sample setter method. It is occasionally syntactically easier to assign samples by reference access.
	 *	\param channel
	 *	\param frame
	 */
	Sample & get_sample( Channel channel, Frame frame );
	
	/** This zero fills the buffer.
	 */
	void clear_buffer();

	/** Raw buffer access. Inderect access is preferred, but for computationally expensive 
	 *	transformations, raw access usually provides an optimization.
	 *	\param channel
	 *	\param frame
	 */
	Sample * get_sample_pointer( Channel channel, Frame frame );

	/** Read-only raw buffer access. Inderect access is preferred, but for computationally expensive 
	 *	transformations, raw access usually provides an optimization.
	 *	\param channel
	 *	\param frame
	 */
	const Sample * get_sample_pointer( Channel channel, Frame frame ) const;

	/** Direct buffer access. Only use this if the data layout doesn't matter or an api requires it.
	 */
	std::vector<Sample> & get_buffer();
	const std::vector<Sample> & get_buffer() const;

	std::vector<Sample>::const_iterator channel_begin( Channel channel ) const;
	std::vector<Sample>::const_iterator channel_end( Channel channel ) const;

	inline size_t get_buffer_pos( Channel, Frame ) const;

#if defined(_WIN32) || defined(WIN32)
	void play() const;
#endif

private://=================================================================================================

	Format format;
	std::vector<Sample> buffer;
};

/** Serialization.
 */
std::ostream & operator<<( std::ostream & os, const AudioBuffer & audio );

} // End namespace flan
