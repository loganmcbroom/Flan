#pragma once

#include <vector>
#include <string>
#include <memory>

#include "flan/defines.h"
#include "flan/Utility/vec2.h"

namespace flan {

/** PVBuffer stores PV data and provides basic buffer access, conversion constants, loading, and saving.
 *	Access to the raw sample buffer is given, but PVBuffer::get_MF and 
 *	PVBuffer::set_MF are preferred when speed is not a factor.
 *	Data is stored in channel -> frame -> bin order, in other words, the entire first channel 
 *	is stored in memory before the second, and within each channel the entire first frame is stored before the second.
 *
 *	What is PV data? PV stands for Phase Vocoder. A phase vocoder transforms audio (or any other signal) into PV data.
 *	The algorithm itself is beyond the scope of this comment, curious readers may find plenty of information on the specifics online.
 *	The output of the algorithm is a two-dimensional set of data, where each data point contains a magnitude and a frequency.
 *	One dimension represents time, the other frequency. These are indexed by "frame" and "bin", respectively.
 *	Both conceptually and literally, each frame of PV data can be seen as an fft of input audio data near that frame's point in time.
 *	The main purpose of the data type in terms of creative audio processing is to allow spectral manipulation as a function of time,
 *	or temporal manipulation as a function of frequency, if you prefer. 
 *	The name "Phase Vocoder" is used only for historic reasons and is a misnomer. 
 */
class PVBuffer
{
public:
	PVBuffer( const PVBuffer & ) = delete;
	PVBuffer( PVBuffer && ) = default;
	PVBuffer& operator=( const PVBuffer & ) = delete;
	PVBuffer& operator=( PVBuffer && ) = default;
	~PVBuffer() = default;

	/** PVBuffer::Format stores all required information outside of the buffer. 
	 *	This is used to easily transfer the way an PVBuffer is stored on disk to a new PVBuffer without copying the buffer.
	 *	In practice, a function transforming a PVBuffer will call PVBuffer::get_format(), modify members as needed, and
	 *	then construct a blank PVBuffer with the Format copy constructor.
	 *	Note that sample rate describes audio frames per second rather than flan frames per second. This is generally easier to work with
	 *	because flan frame rates can be non-integer values.
	 */
	struct Format 
		{ 
		Channel num_channels = 0;
		Frame num_frames = 0;
		Bin num_bins = 0;
		FrameRate sample_rate = 48000; 
		FrameRate analysis_rate = 48000 / 128;
		Frame window_size = 0;
		//window type?
		};

	/** Default constructor
	 */
	PVBuffer();

	/** Format copy constructor, constructs a PVBuffer with the given format and an uninitialized buffer.
	 *	\param format Format to use in construction.
	 */
	PVBuffer( const Format & format );

	/** Load constructor, calls PVBuffer::load.
	 *	\param filepath A PV file to load.
	 */
	PVBuffer( const std::string & filename );

	/** Returns a deep copy of the PVBuffer.
	 */
	PVBuffer copy() const;

	/** Return true if the PV is in a state which cannot be processed. This includes a 0-channel buffer,
	 *	a 0-frame buffer, a 0-bin buffer, or a buffer with a 0 sample rate.
	 */
	bool is_null() const;

	/// @brief Scans all samples for nan or inf.
	/// @return Returns true if any sample is nan or inf.
	bool is_nan_or_inf() const;

	//======================================================
	//	I/O
	//======================================================

	/** File loading. This utilizes the flan specific RIFF data type, PV, which is defined here. PV files should use the extension flan.
	 *	Data should be saved in little-endian format. 
	 *
	 *	Chunk one is the RIFF chunk. 
	 *		Bytes 0-3 is "RIFF". 
	 *		Bytes 4-7 is 4 (uint32_t), the size of the RIFF chunk (12) minus the data up to and including this data (8).
	 *		Bytes 8-11 is "PV".
	 *
	 *	Chunk two is the PV format chunk.
	 *		Bytes 0-3 is "fmt " (note the space).
	 *		Bytes 4-7 is 26 (uint32_t), the size of the format chunk (34) minus the data up to and including this data (8).
	 *		Bytes 8-9 is 1 (uint16_t), this indicates default formatting. This is unused and reserved for future formatting changes.
	 *		Bytes 10-11 is the number of channels (uint16_t).
	 *		Bytes 12-15 is the number of frames (uint32_t).
	 *		Bytes 16-19 is the number of bins per frame (uint32_t).
	 *		Bytes 20-23 is the audio sample rate used to create the PV data. This is used because the PV frame rate can be a non-integer value.
	 *		Bytes 24-27 is the hop size used in the phase vocoder (uint32_t).
	 *		Bytes 28-31 is the number of bytes used to store buffer information (uint32_t). Each MF pair stores twice this number of bytes.
	 *		Bytes 32-33 is a phase vocoder window function id. Currently only 1 is defined and represents a hann window.
	 *	
	 *	Chunk three is the data chunk.
	 *		Bytes 0-3 is "data".
	 *		Bytes 4-7 is the size of the data chunk minus the data up to and including this data (8).
	 *		The remainder of the bytes are the flan data stored in channel -> frame -> bin order. 
	 *		Each piece of PV data should be stored in magnitude, frequency order, using signed integers (as WAVE does).
	 *		Note that magnitudes will likely need to be normalized before storage.
	 *		Frequency data should be scaled so the max 24bit signed int value maps to the audio sample rate corresponding to the PV saved.
	 *
	 *	\param filepath A PV file to load.
	 */
	bool load( const std::string & filename );

	/** File saving. See PVBuffer::load for format information.
	 *	\param filepath A PV file to load.
	 */
	bool save( const std::string & filename ) const;

	/** Prints format data.
	 */
	void print_summary() const;

	//======================================================
	//	Getters
	//======================================================

	/** Main MF access method.
	 *	\param channel
	 *	\param frame
	 *	\param bin
	 */
	MF get_MF( Channel channel, Frame frame, Bin bin ) const;

	/** Returns a single frame PVBuffer with data from the given input frame.
	 * \param frame
	 */
	PVBuffer get_frame( Frame frame ) const;

	/** Returns the format being used.
	 */
	Format get_format() const;

	/** Returns the number of channels.
	 */
	Channel get_num_channels() const;

	/** Returns the number of frames.
	 */
	Frame get_num_frames() const;

	/** Returns the number of bins.
	 */
	Bin get_num_bins() const;

	/** Returns the frames per second of the audio from which this PVBuffer came.
	 *	If the frames per second of the PV data is needed, use PVBuffer::get_analysis_rate.
	 */
	FrameRate get_sample_rate() const;

	/** Returns the spectral frames per second.
	 *	If the frames per second of the auto data is needed, use PVBuffer::get_sample_rate.
	 */
	FrameRate get_analysis_rate() const;

	/** Returns the hop size used in the phase vocoder to create this PVBuffer.
	 */
	Frame get_hop_size() const;

	/** Returns the size of the transform used in the phase vocoder to create this PVBuffer.
	 */
	Frame get_dft_size() const;

	/** Returns the size of the input data window used in the phase vocoder to create this PVBuffer.
	 */
	Frame get_window_size() const;

	/** Returns length in seconds.
	 */
	Second get_length() const;

	/** Returns height in hz.
	 */
	Frequency get_height() const;

	/** Returns the magnitude of whichever MF has the largest magnitude.
	 */
	Magnitude get_max_partial_magnitude() const;

	/** Returns the magnitude of whichever MF has the largest magnitude in the given range.
	 * \param start_frame The start of the time range.
	 * \param end_frame The end of the time range. Passing 0 will use the maximum time.
	 * \param start_bin The start of the frequency range.
	 * \param end_bin The end of the frequency range. Passing 0 will use the maximum frequency.
	 */
	Magnitude get_max_partial_magnitude( uint32_t start_frame, uint32_t end_frame = 0, 
		uint32_t start_bin = 0, uint32_t end_bin = 0 ) const;

	/** Returns a unit fraction for converting seconds to frames
	 */
	fFrame time_to_frame( Second ) const;

	/** Returns a unit fraction for converting frames to seconds.
	 */
	Second frame_to_time( fFrame ) const;

	/** Returns a unit fraction for converting frequencies to bins.
	 */
	fBin frequency_to_bin( Frequency ) const;

	/** Returns a unit fraction for converting bins to frequencies.
	 */
	Frequency bin_to_frequency( fBin ) const;

	/** Returns how far the frequency recorded at the given position is from the bin frequency.
	 */
	Frequency get_frequency_offset( Channel c, Frame f, Bin b ) const;

	/** Returns the input channel clamped to the buffer bounds
	 */
	Channel bound_channel( Channel c ) const;

	/** Returns the input frame clamped to the buffer bounds
	 */
	Frame bound_frame( Frame c ) const;

	/** Returns the input bin clamped to the buffer bounds
	 */
	Bin bound_bin( Bin c ) const;

	//======================================================
	//	Setters
	//======================================================

	/** Main MF setter method.
	 *	\param channel
	 *	\param frame
	 *	\param bin
	 *	\param mf The value to write.
	 */
	void set_MF( Channel channel, Frame frame, Bin bin, MF mf );

	/** Additional MF setter method. It is occasionally syntactically easier to assign MFs by reference access.
	 *	\param channel
	 *	\param frame
	 *	\param bin
	 */
	MF & get_MF( Channel channel, Frame frame, Bin bin );

	/** This zero fills the buffer.
	 */
	void clear_buffer();

	/** Raw buffer access. Inderect access is preferred, but for computationally expensive 
	 *	transformations, raw access usually provides an optimization.
	 *	\param channel
	 *	\param frame
	 *	\param bin
	 */
	MF * get_MF_pointer( Channel channel, Frame frame, Bin bin );

	/** Read-only raw buffer access. Inderect access is preferred, but for computationally expensive 
	 *	transformations, raw access usually provides an optimization.
	 *	\param channel
	 *	\param frame
	 *	\param bin
	 */
	const MF * get_MF_pointer( Channel channel, Frame frame, Bin bin  ) const;

	std::vector<MF> & get_buffer();
	const std::vector<MF> & get_buffer() const;

	std::vector<MF>::iterator channel_begin( Channel channel );
	std::vector<MF>::iterator channel_end( Channel channel );
	std::vector<MF>::const_iterator channel_begin( Channel channel ) const;
	std::vector<MF>::const_iterator channel_end( Channel channel ) const;

	/** This converts from a channel, frame, and bin, to the appropriate position in the buffer.
	 */
	size_t get_buffer_pos( Channel, Frame, Bin ) const;

private: //=================================================================================================
	
	PVBuffer::Format format;
	std::vector<MF> buffer;
};

/** Serialization.
 */
std::ostream & operator<<( std::ostream & os, const PVBuffer & flan );

} // End namespace flan

