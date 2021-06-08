#pragma once

#include <vector>
#include <string>
#include <memory>

#include "defines.h"

namespace flan {

/** PVOCBuffer stores PVOC data and provides basic buffer access, conversion constants, loading, and saving.
 *	Access to the raw sample buffer is given, but PVOCBuffer::getMF and 
 *	PVOCBuffer::setMF are preferred when speed is not a factor.
 *	Data is stored in channel -> frame -> bin order, in other words, the entire first channel 
 *	is stored in memory before the second, and within each channel the entire first frame is stored before the second.
 *
 *	What is PVOC data? PVOC stands for Phase Vocoder. A phase vocoder transforms audio (or any other signal) into PVOC data.
 *	The algorithm itself is beyond the scope of this comment, curious readers may find plenty of information on the specifics online.
 *	The output of the algorithm is a two-dimensional set of data, where each data point contains a magnitude and a frequency.
 *	One dimension represents time, the other frequency. These are indexed by "frame" and "bin", respectively.
 *	Both conceptually and literally, each frame of PVOC data can be seen as an fft of input audio data near that frame's point in time.
 *	The main purpose of the data type in terms of creative audio processing is to allow spectral manipulation as a function of time,
 *	or temporal manipulation as a function of frequency, if you prefer. 
 *	The name "Phase Vocoder" is used only for historic reasons and is a misnomer. 
 */
class PVOCBuffer
{
public:

	/** This is the data type stored at each buffer position. It represents a frequency, and the magnitude of that frequency.
	*/
	struct MF 
		{ 
		Magnitude m;
		Frequency f; 
		};

	/** PVOCBuffer::Format stores all required information outside of the buffer. 
	 *	This is used to easily transfer the way an PVOCBuffer is stored on disk to a new PVOCBuffer without copying the buffer.
	 *	In practice, a function transforming a PVOCBuffer will call PVOCBuffer::getFormat(), modify members as needed, and
	 *	then construct a blank PVOCBuffer with the Format copy constructor.
	 *	Note that sample rate describes audio frames per second rather than flan frames per second. This is generally easier to work with
	 *	because flan frame rates can be non-integer values.
	 */
	struct Format 
		{ 
		Channel numChannels = 0;
		Frame numFrames = 0;
		Bin numBins = 0;
		uint32_t sampleRate = 48000; 
		Frame hopSize = 0;
		Frame windowSize = 0;
		//window type?
		};

	/** Default constructor
	 */
	PVOCBuffer();

	/** Format copy constructor, constructs a PVOCBuffer with the given format and an uninitialized buffer.
	 *	\param format Format to use in construction.
	 */
	PVOCBuffer( const Format & );

	/** Load constructor, calls PVOCBuffer::load.
	 *	\param filePath A PVOC file to load.
	 */
	PVOCBuffer( const std::string & filename );

	/** Returns a deep copy of the PVOCBuffer.
	 */
	PVOCBuffer deepCopy() const;

	/** Return true if the PVOC is in a state which cannot be processed. This includes a 0-channel buffer,
	 *	a 0-frame buffer, a 0-bin buffer, or a buffer with a 0 sample rate.
	 */
	bool isNull() const;

	//======================================================
	//	I/O
	//======================================================

	/** File loading. This utilizes the flan specific RIFF data type, PVOC, which is defined here. PVOC files should use the extension flan.
	 *	Data should be saved in little-endian format. 
	 *
	 *	Chunk one is the RIFF chunk. 
	 *		Bytes 0-3 is "RIFF". 
	 *		Bytes 4-7 is 4 (uint32_t), the size of the RIFF chunk (12) minus the data up to and including this data (8).
	 *		Bytes 8-11 is "PVOC".
	 *
	 *	Chunk two is the PVOC format chunk.
	 *		Bytes 0-3 is "fmt " (note the space).
	 *		Bytes 4-7 is 26 (uint32_t), the size of the format chunk (34) minus the data up to and including this data (8).
	 *		Bytes 8-9 is 1 (uint16_t), this indicates default formatting. This is unused and reserved for future formatting changes.
	 *		Bytes 10-11 is the number of channels (uint16_t).
	 *		Bytes 12-15 is the number of frames (uint32_t).
	 *		Bytes 16-19 is the number of bins per frame (uint32_t).
	 *		Bytes 20-23 is the audio sample rate used to create the PVOC data. This is used because the PVOC frame rate can be a non-integer value.
	 *		Bytes 24-27 is the hop size used in the phase vocoder (uint32_t).
	 *		Bytes 28-31 is the number of bytes used to store buffer information (uint32_t). Each MF pair stores twice this number of bytes.
	 *		Bytes 32-33 is a phase vocoder window function id. Currently only 1 is defined and represents a Hann window.
	 *	
	 *	Chunk three is the data chunk.
	 *		Bytes 0-3 is "data".
	 *		Bytes 4-7 is the size of the data chunk minus the data up to and including this data (8).
	 *		The remainder of the bytes are the flan data stored in channel -> frame -> bin order. 
	 *		Each piece of PVOC data should be stored in magnitude, frequency order, using signed integers (as WAVE does).
	 *		Note that magnitudes will likely need to be normalized before storage.
	 *		Frequency data should be scaled so the max 24bit signed int value maps to the audio sample rate corresponding to the PVOC saved.
	 *
	 *	\param filePath A PVOC file to load.
	 */
	bool load( const std::string & filename );

	/** File saving. See PVOCBuffer::load for format information.
	 *	\param filePath A PVOC file to load.
	 */
	bool save( const std::string & filename ) const;

	/** Prints format data.
	 */
	void printSummary() const;

	//======================================================
	//	Getters
	//======================================================

	/** Main MF access method.
	 *	\param channel
	 *	\param frame
	 *	\param bin
	 */
	MF getMF( Channel channel, Frame frame, Bin bin ) const;

	/** Returns a single frame PVOCBuffer with data from the given input frame.
	 * \param frame
	 */
	PVOCBuffer getFrame( Frame frame ) const;

	/** Returns the format being used.
	 */
	Format getFormat() const;

	/** Returns the number of channels.
	 */
	Channel getNumChannels() const;

	/** Returns the number of frames.
	 */
	Frame getNumFrames() const;

	/** Returns the number of bins.
	 */
	Bin getNumBins() const;

	/** Returns the frames per second of the audio from which this PVOCBuffer came.
	 *	If the frames per second of the PVOC data is needed, use PVOCBuffer::timeToFrame.
	 */
	uint32_t getSampleRate() const;

	/** Returns the hop size used in the phase vocoder to create this PVOCBuffer.
	 */
	Frame getHopSize() const;

	/** Returns the size of the transform used in the phase vocoder to create this PVOCBuffer.
	 */
	Frame getDFTSize() const;

	/** Returns the size of the input data window used in the phase vocoder to create this PVOCBuffer.
	 */
	Frame getWindowSize() const;

	/** Returns length in seconds.
	 */
	Time getLength() const;

	/** Returns height in hz.
	 */
	Frequency getHeight() const;

	/** Returns the magnitude of whichever MF has the largest magnitude.
	 */
	Magnitude getMaxPartialMagnitude() const;

	/** Returns the magnitude of whichever MF has the largest magnitude in the given range.
	 * \param startFrame The start of the time range.
	 * \param endFrame The end of the time range. Passing 0 will use the maximum time.
	 * \param startBin The start of the frequency range.
	 * \param endBin The end of the frequency range. Passing 0 will use the maximum frequency.
	 */
	Magnitude getMaxPartialMagnitude( uint32_t startFrame, uint32_t endFrame = 0, 
		uint32_t startBin = 0, uint32_t endBin = 0 ) const;

	/** Returns a unit fraction for converting seconds to frames
	 */
	float timeToFrame() const;

	/** Returns a unit fraction for converting frames to seconds.
	 */
	float frameToTime() const;

	/** Returns a unit fraction for converting frequencies to bins.
	 */
	float frequencyToBin() const;

	/** Returns a unit fraction for converting bins to frequencies.
	 */
	float binToFrequency() const;

	/** Returns how far the frequency recorded at the given position is from the bin frequency.
	 */
	float getFrequencyOffset( Channel c, Frame f, Bin b ) const;

	/** Returns the input channel clamped to the buffer bounds
	 */
	Channel boundChannel( Channel c ) const;

	/** Returns the input frame clamped to the buffer bounds
	 */
	Frame boundFrame( Frame c ) const;

	/** Returns the input bin clamped to the buffer bounds
	 */
	Bin boundBin( Bin c ) const;

	//======================================================
	//	Setters
	//======================================================

	/** Main MF setter method.
	 *	\param channel
	 *	\param frame
	 *	\param bin
	 *	\param mf The value to write.
	 */
	void setMF( Channel channel, Frame frame, Bin bin, MF mf );

	/** Additional MF setter method. It is occasionally syntactically easier to assign MFs by reference access.
	 *	\param channel
	 *	\param frame
	 *	\param bin
	 */
	MF & getMF( Channel channel, Frame frame, Bin bin );

	/** This zero fills the buffer.
	 */
	void clearBuffer();

	/** Raw buffer access. Inderect access is preferred, but for computationally expensive 
	 *	transformations, raw access usually provides an optimization.
	 *	\param channel
	 *	\param frame
	 *	\param bin
	 */
	MF * getMFPointer( Channel channel, Frame frame, Bin bin );

	/** Read-only raw buffer access. Inderect access is preferred, but for computationally expensive 
	 *	transformations, raw access usually provides an optimization.
	 *	\param channel
	 *	\param frame
	 *	\param bin
	 */
	const MF * getMFPointer( Channel channel, Frame frame, Bin bin  ) const;

	std::vector<MF>::iterator begin();
	std::vector<MF>::iterator end();
	std::vector<MF>::const_iterator begin() const;
	std::vector<MF>::const_iterator end() const;

	/** This converts from a channel, frame, and bin, to the appropriate position in the buffer.
	 */
	size_t getBufferPos( Channel, Frame, Bin ) const;

private: //=================================================================================================
	
	Format format;
	std::shared_ptr<std::vector< MF >> buffer;
};

/** Serialization.
 */
std::ostream & operator<<( std::ostream & os, const PVOCBuffer & flan );

} // End namespace flan

