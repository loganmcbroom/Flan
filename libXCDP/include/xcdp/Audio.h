#pragma once

#include <functional>

#include "AudioBuffer.h"

namespace xcdp {

class PVOC;
class Spectrum;
struct Func1x1;

//! Audio Algorithms
/** This is a stateless wrapper of xcdp::AudioBuffer which contains all relevant algorithms and utilities available
 *	in xcdp for manipulating Audio buffers.
 *
 *	All methods are const and have no side effects outside of Audio::graph, which saves a bmp file.
 */
class Audio : public AudioBuffer
{
public:
	typedef const std::vector<Audio> & Vec;

	/** Iterative methods can take these, usually passing them the input Audio and the current iteration number.
	 */
	typedef std::function< Audio ( const Audio &, size_t ) > Mod;
	
	/** Constructs 0 size Audio with default AudioBuffer::Format. 
	 */
	Audio();

	/** Constructs an Audio with the Format other but with an initialized buffer.
	 *	\param other Format for the constructed PVOC
	 */
	Audio( const AudioBuffer::Format & other );

	/** Constructs an Audio from the file at filename.
	 *	Supports any file extensions supported by libsndfile. Wave is suggested, and mp3 is unsupported.
	 *	\param filename Filename to load
	 */
	Audio( const std::string & filename );

	/** Helper for simple mixing of two Audio, see Audio::mix.
	 *	\param other Audio to mix with this
	 */
	Audio operator+( const Audio & other ) const;

	/** Helper for phase inversion, see Audio::invertPhase
	 */
	Audio operator-() const;

	/** Find the total energy in the input, in other words, the sum of the square of each sample. 
	 */
	float getTotalEnergy() const;

	/** Creates and saves a bmp of the input waveform. Each channel is graphed and the images are stacked vertically.
	 *	\param filename File path at which to save (should end in .bmp)
	 *	\param width Output bmp width
	 *	\param height The height of each channel in the bmp.
	 */
	Audio graph( const std::string & filename, size_t width = 2048, size_t height = 512 ) const;

	/** This is a phase vocoder, the main transform from Audio to PVOC data. 
	 *	See PVOCBuffer for information on the meaning and structure of the PVOC data type.
	 *	This algorithm is implemented using OpenCL by default. 
	 *	If you are building without OpenCL, this will call a cpu implementation.
	 *
	 *	\param frameSize This is the number of samples taken for each fft. 
	 *		Defined behaviour is only guaranteed for power-of-two inputs, but it may work for other input sizes.
	 *		The number of bins in the output PVOC is given by frameSize / 2 + 1.
	 *	\param overlaps Taking frameSize / overlaps gives the forward jump in the input Audio per fft step.
	 *		This can be thought of as a time-oversampling factor. 
	 *		The size (total data) of the output PVOC is the product of the input audio size, two, and overlaps.
	 */
	PVOC convertToPVOC( const size_t frameSize = 2048, const size_t overlaps = 16 ) const;

	//Spectrum convertToSpectrum() const;

	/** AudioBuffer normally stores "Left" and "Right" data in channel 0 and 1 respectively.
	 *	The result of this transform stores "Mid" and "Side" data in those channels instead.
	 *	This allows some simple effects such as widening, but also allows mid or side information
	 *	to be processed indepentantly of the other.
	 */
	Audio convertToMidSide() const;

	/** This converts back to Left/Right from Mid/Size. See Audio::convertToMidSide. 
	 */
	Audio convertToLeftRight() const;

	/** This converts the input to a two channel Audio. The conversion is channel-count dependant.
	 *	The currently accepted channel counts are: 1 and 2.
	 */
	Audio convertToStereo() const;

	/** This converts an Audio with any number of channels to a single channel Audio.
	 *	Each output sample is the sum of the channels at that time over the square root of the number of channels.
	 */
	Audio convertToMono() const;

	/** This phase inverts the input. Every output sample is assigned the negative of the input sample at that time.
	 */
	Audio invertPhase() const;

	/** Basic volume scaling, Output( t ) = input( t ) * volumeLevel( t ).
	 *	\param gain Describes how volume should be scaled as a function of time.
	 */
	Audio modifyVolume( Func1x1 gain ) const;

	/** Basic volume setting. This normalizes (divides every sample by the largest sample value), and then scales by level.
	 *	\param level Describes how volume should be scaled as a function of time.
	 */
	Audio setVolume( Func1x1 level ) const;

	/** This applies the shaper as a function to each sample in the input.
	 *	\param shaper Each sample in the input is passed through this. 
	 *		Samples will be values on [-1,1] under normal circumstances.
	 *		For example, the function y = x would have no effect as a shaper.
	 */
	Audio waveshape( Func1x1 shaper ) const;

	/** This spatially repositions the input. Panning is channel-count dependant.
	 *	\param lrAmount This should return a value on [-1,1], representing movement from left to right.
	 */
	Audio pan( Func1x1 lrAmount ) const;

	/** This redistributes energy between the mid and side signals
	 *	\param widenAmount This should return a value on [-1,1], representing movement from mid to side.
	 */
	Audio widen( Func1x1 widenAmount ) const;

	/** This repeats the input n times. Each iteration is first fed into mod along with which iteration is being sent.
	 *	\param n The number of desired iterations
	 *	\param mod A function taking an Audio and the iteration index.
	 *	\param fbIterate When this flag is set, the previous iteration is sent to mod, rather than the original Audio.
	 */
	Audio iterate( size_t n, Audio::Mod mod = nullptr, bool fbIterate = false ) const;

	/** This reverses the Audio.
	 */
	Audio reverse() const;

	/** This returns a peice of the original Audio that lied between startTime and endTime.
	 *	\param startTime Start time for cut.
	 *	\param endTime End time for cut.
	 */
	Audio cut( float startTime, float endTime ) const;

	/** This repitches the input.
	 *	\param factor The output of this represents a scaling factor of pitch as a function of time.
	 *	\param granularity This represents the frequency at which factor is sampled, in number of samples.
			Rapidly changing factors should use a small granularity.
	 *	\param qual This can take the value 0 for a default quality sinc resampling, 1 for quick and dirty (lerp and simple filter), 
	 *		or 2 for bad quality.
	 */
	Audio repitch( Func1x1 factor, size_t granularity = 64, size_t qual = 0 ) const;

	/** If each Func1x1 in ir is constant this is equivalent to normal convolution between the input Audio and ir.
	 *	\param ir An impulse response, this is treated like an Audio whos samples can change over time. 
	 */
	Audio convolve( const std::vector<Func1x1> & ir ) const;

	/** This is normal convolution between Audio.
	 *	\param ir An impulse response.
	 */
	Audio convolve( const Audio & ir ) const;

	/** Basic delay effect. Due to the non-real-time nature of xcdp, this can only produce a finite output, and 
	 *	thus a set number of delays must be specified.
	 *	\param delayTime The amount of time between delays.
	 *	\param numDelays The number of delays, 1 does nothing.
	 *	\param decay The nth decay has its volume scaled by decay^n.
	 *	\param mod This allows user-defined processing of each delay.
	 *	\param fbIterate When this flag is set, the previous delay is sent to mod, rather than the original Audio.
	 */
	Audio delay( float delayTime, size_t numDelays, float decay = .5, Audio::Mod mod = nullptr, bool fbIterate = true ) const;

	/** This adds a cosine based fade to the ends of the input Audio.
	 *	\param fadeTime The length of the fades.
	 *	\param start This flag enables fading in the start of the input Audio.
	 *	\param end This flag enables fading out the end of the input Audio.
	 */
	Audio fades( float fadeTime = .05, bool start = true, bool end = true ) const;

	/** A simple low pass filter. This method will likely be overhauled with the planned filter methods.
	 *	\param cutoff The frequency at which attenuation should begin as a function of time.
	 *	\param taps This in a sense describes the "accuracy" of the filter. 
			Increasing it will produce a sharper cutoff and increase processing time.
	 */
	Audio lowPass( Func1x1 cutoff, size_t taps = 64 ) const;

	// Unfinished
	//Audio freeze( real_t freezeTime, real_t freezeLength, real_t delay, real_t rand = 0, Audio::Mod mod = nullptr ) const;
	// drunk walk / scramble
	// extend loop
	// extend repititions?

	//========================================================
	// Multi-In Procs
	//========================================================

	/** The main Audio mixing method.
	 *	\param ins The Audio to mix.
	 *	\param balances The volume scaling to apply to each input Audio. 
	 *		The default input will balance all inputs to 1/n, where n is the number of inputs.
	 *	\param startTimes This gives the start time for each input Audio. The default will start all inputs at time 0.
	 */
	static Audio mix( Audio::Vec ins, 
		std::vector<Func1x1> balances = std::vector<Func1x1>(),
		std::vector<float> startTimes = std::vector<float>() );

	/** This joins all input Audio tip to tail in the order they were passed.
	 *	\param ins The Audio to mix.
	 */
	static Audio join( Audio::Vec ins );
};

} // End namespace xcdp
