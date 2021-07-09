#pragma once

#include <functional>

#include "AudioBuffer.h"
#include "flan/Utility/Interpolator.h"
#include "flan/Utility/Interval.h"
#include "Function.h"

namespace std { template<class T> class complex; };

namespace flan {

class PVOC;
class Spectrum;
class FFTHelper;
class Graph;

//! Audio Algorithms
/** This is a wrapper of flan::AudioBuffer which contains all relevant algorithms and utilities available
 *	in flan for manipulating Audio buffers.
 *
 *	All Audio methods are const and have no side effects.
 */
class Audio : public AudioBuffer
{
public:
	/** Iterative methods utilize these, usually passing them the input Audio and time.
	 */
	using Mod = std::function< Audio ( const Audio &, Time ) >;
	
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

	/** Audio is a stateless wrapper, so it can be constructed from a buffer
	 */
	Audio( const AudioBuffer & );

	/** Helper for simple mixing of two Audio, see Audio::mix.
	 *	\param other Audio to mix with this
	 */
	Audio operator+( const Audio & other ) const;

	/** Helper for phase inversion, see Audio::invertPhase
	 */
	Audio operator-() const;

	/** Helper for finding Audio differences. Returns *this + -other.
	 */
	Audio operator-( const Audio & other ) const;

	/** Find the total energy in the input, in other words, the sum of the square of each sample. 
	 */
	float getTotalEnergy( flan_CANCEL_ARG ) const;

	/** Access the Audio buffer as a continuous function by interpolating between sample values.
	 *	 Note, this interpolation is not "smooth" in the way a resampler interpolates, it is a basic
	 *	 interpolation considering only two points at a time.
	 *	\param channel The channel to access.
	 *	\param time The time to access.
	 *	\param interp How the interpolation should be computed.
	 */
	Sample getSampleInterpolated( uint32_t channel, float frame, Interpolator interp = Interpolators::linear ) const;



	//============================================================================================================================================================
	// Conversions
	//============================================================================================================================================================

	/** Converts the Audio to a waveform bmp.
	 *  \param tStart The graph start time.
	 *  \param tEnd The graph end time. Passing -1 will graph up to the end of the data.
	 *  \param width The bmp width.
	 *  \param height The bmp height.
	 *	\param timelineScale Size of tick marks. Zero for none.
	 */
	Graph convertToGraph( Interval I = Interval( 0, -1 ), Pixel width = -1, Pixel height = -1, float timelineScale = 0, flan_CANCEL_ARG ) const;

	/** Converts the Autio to a waveform bmp and saves it at filename.
	 *	\param filename A filename at which to save the generated bmp.
	 *  \param I The interval of time to graph.
	 *  \param width The bmp width.
	 *  \param height The bmp height.
	 */
	Audio saveToBMP( const std::string & filename, Interval I = Interval( 0, -1 ), Pixel width = -1, Pixel height = -1, flan_CANCEL_ARG ) const;

	/** This is a phase vocoder, the main transform from Audio to PVOC data. 
	 *	See PVOCBuffer for information on the meaning and structure of the PVOC data type.
	 *	This algorithm is implemented using OpenCL by default. 
	 *	If you are building without OpenCL, this will call a cpu implementation.
	 *
	 *	\param windowSize This is the number of Frames copied into each fft input. 
	 *		Defined behaviour is only guaranteed for power-of-two inputs, but it may work for other input sizes.
	 *		The number of bins in the output PVOC is given by frameSize / 2 + 1.
	 *	\param hop This is the forward jump in the input Audio per fft step.
	 *  \param fftSize The size of the fft. 
	 *  \param fftwInBuffer When fftw is used this allows allocating the fftw input buffer used in the transform
	 *		on the main thread, as fftw buffer allocations aren't thread safe.
	 *  \param fftwOutBuffer When fftw is used this allows allocating the fftw output buffer used in the transform
	 *		on the main thread, as fftw buffer allocations aren't thread safe.
	 *	 \param plan When fftw is used this allows allocating the fftw plan used in the transform
	 *		on the main thread, as fftw planning isn't thread safe.
	 */
	PVOC convertToPVOC( Frame windowSize = 2048, Frame hop = 128, Frame fftSize = 4096, std::shared_ptr<FFTHelper> fft = nullptr, flan_CANCEL_ARG ) const;

	/** See Audio::convertToPVOC. This explicitly calls a cpu-based implementation.
	 */
	PVOC convertToPVOC_cpu( Frame windowSize = 2048, Frame hop = 128, Frame fftSize = 4096, std::shared_ptr<FFTHelper> fft = nullptr, flan_CANCEL_ARG ) const;

	/** This extracts a function approximating volume over time. The extracted function 
	 *	is generated by linear approximations. Only the first channel is processed.
	 *
	 *  \param granularity This parameter decides how often (in time) the current average gain is sampled. 
	 *		A granularity of 1 will return a function exactly equal to the input Audio.
	 *		Smaller granularities will generate more costly functions.
	 */
	Func1x1 convertToFunction( Time granularity = .001f, flan_CANCEL_ARG ) const;

	/** AudioBuffer normally stores "Left" and "Right" data in channel 0 and 1 respectively.
	 *	The result of this transform stores "Mid" and "Side" data in those channels instead.
	 *	This allows some simple effects such as widening, but also allows mid and side information
	 *	to be processed independantly.
	 */
	Audio convertToMidSide( flan_CANCEL_ARG ) const;

	/** This converts back to Left/Right from Mid/Side. See Audio::convertToMidSide. 
	 */
	Audio convertToLeftRight( flan_CANCEL_ARG ) const;

	/** This converts the input to a two channel Audio. The conversion is channel-count dependant.
	 *	The currently accepted channel counts are: 1 and 2.
	 */
	Audio convertToStereo( flan_CANCEL_ARG ) const;

	/** This converts an Audio with any number of channels to a single channel Audio.
	 *	Each output sample is the sum of the channels at that time over the square root of the number of channels.
	 */
	Audio convertToMono( flan_CANCEL_ARG ) const;



	//============================================================================================================================================================
	// Information
	//============================================================================================================================================================

	/** Try to find the wavelength of the input over time. Selects to minimize differences per repitition.
	 *	\param channel The channel to process.
	 *	\param start Frame to analyze.
	 *	\param absoluteCutoff Minimum the local wave difference minimum must achieve to be considered a potential wavelength.
	 *	\param windowSize Input window size.
	 *	\param fft External fft handler. IMPORTANT: if using this, the fft size needs to be found with autocorrelationFFTSize( windowSize ) in DSPUtility.
	 */
	float getLocalWavelength( Channel channel, Frame start, Frame windowSize = 2048, std::shared_ptr<FFTHelper> fft = nullptr, flan_CANCEL_ARG ) const;

	/** Try to find the wavelength of the input over time. This is estimated using parabolic interpolation, so it may not be an integer.
	 *	\param channel The channel to process.
	 *	\param start Start frame.
	 *	\param end End frame.
	 *	\param absoluteCutoff Minimum the local wave difference minimum must achieve to be considered a potential wavelength.
	 *	\param windowSize The amount of input data used in each analysis.
	 *	\param hop The hop size in frames per anlysis.
	 *	\param fft The fft object. See flan::FFTHelper.
	 */
	std::vector<float> getLocalWavelengths( Channel channel, Frame start = 0, Frame end = -1, Frame windowSize = 2048, Frame hop = 128, std::shared_ptr<FFTHelper> fft = nullptr, flan_CANCEL_ARG ) const;

	/** Find the average of local wavelengths over time.
	 *	\param localWavelengths Wavelengths to average. A -1 wavelength indicates no detectable wavelength. 
	 *	\param minActiveRatio This proportion of the inputs must have a detectable wavelength to continue computation. Otherwise, -1 is returned.
	 *	\param maxLengthSigma If the standard deviation of the wavelengths is above this threshold, -1 is returned. Inputting -1 will disable sigma thresholding.
	 */
	float getAverageWavelength( const std::vector<float> & localWavelengths, float minActiveRatio = 0, float maxLengthSigma = -1 ) const; // Same as getAverageWavelengths

	/** Find the average of local wavelengths over time. This is a utility for calling getLocalWavelength and passing it to the other getAverageWavelength.
	 *	\param channel The channel to process.
	 *	\param minActiveRatio This proportion of the inputs must have a detectable wavelength to continue computation. Otherwise, -1 is returned.
	 *	\param maxLengthSigma If the standard deviation of the wavelengths is above this threshold, -1 is returned. Inputting -1 will disable sigma thresholding.
	 *	\param start Start frame.
	 *	\param end End frame.
	 *	\param windowSize The amount of input data used in each analysis.
	 *	\param hop The hop size in frames per anlysis.
	 *	\param fft The fft object. See flan::FFTHelper.
	 */
	float getAverageWavelength( Channel, float minActiveRatio = 0, float maxLengthSigma = -1, Frame start = 0, Frame end = -1, 
		Frame windowSize = 2048, Frame hop = 128, std::shared_ptr<FFTHelper> fft = nullptr, flan_CANCEL_ARG ) const; // Averages getLocalWavelengths
	

	/** Try to find the wavelength of the input over time. Selects highest to minimize octave error.
	 *	\param channel The channel to process.
	 *	\param start Frame to analyze.
	 *	\param windowSize Input window size.
	 *	\param fft External fft handler. IMPORTANT: if using this, the fft size needs to be found with autocorrelationFFTSize( windowSize ) in DSPUtility.
	 */
	Frequency getLocalFrequency( Channel channel, Frame start = 0, Frame windowSize = 2048, std::shared_ptr<FFTHelper> fft = nullptr, flan_CANCEL_ARG ) const;

	/** Utility for calling getLocalFrequency at a set frame hop interval from start to end.
	 */
	std::vector<Frequency> getLocalFrequencies( Channel channel, Frame start = 0, Frame end = -1, Frame windowSize = 2048, 
		Frame hop = 128, std::shared_ptr<FFTHelper> fft = nullptr, flan_CANCEL_ARG ) const;

	//============================================================================================================================================================
	// Processing
	//============================================================================================================================================================

	/** This phase inverts the input. Every output sample is assigned the negative of the input sample at that time.
	 */
	Audio invertPhase( flan_CANCEL_ARG ) const;

	/** Volume scaling, Output( t ) = input( t ) * volumeLevel( t ).
	 *	\param gain Describes how volume should be scaled as a function of time.
	 */
	Audio modifyVolume( Func1x1 gain, flan_CANCEL_ARG ) const;

	/** Volume setting. This normalizes (divides every sample by the largest sample value), and then scales by level.
	 *	\param level Describes how volume should be scaled as a function of time.
	 */
	Audio setVolume( Func1x1 level, flan_CANCEL_ARG ) const;

	/** Time shifting.
	 *	\param shift How much to shift the input.
	 */
	Audio shift( Time shift, flan_CANCEL_ARG ) const;

	/** This applies the shaper as a function to each sample in the input.
	 *	\param shaper Each sample in the input is passed through this. 
	 *		Samples will be values on [-1,1] under normal circumstances.
	 *		For example, the function y = x would have no effect as a shaper.
	 */
	Audio waveshape( Func1x1 shaper, flan_CANCEL_ARG ) const;

	/** This spatially repositions the input. Panning is channel-count dependant.
	 *	\param lrAmount This should return a value on [-1,1], representing movement from left to right.
	 */
	Audio pan( Func1x1 lrAmount, flan_CANCEL_ARG ) const;

	/** This redistributes energy between the mid and side signals
	 *	\param widenAmount This should return a value on [-1,1], representing movement from mid to side.
	 */
	Audio widen( Func1x1 widenAmount, flan_CANCEL_ARG ) const;

	/** This reverses the Audio.
	 */
	Audio reverse( flan_CANCEL_ARG ) const;

	/** This returns a peice of the original Audio that lied between startTime and endTime.
	 *	\param start Start time for cut.
	 *	\param end End time for cut.
	 *	\param startFade The cut audio start is faded for this amount of time.
	 *	\param endFade The cut audio end is faded for this amount of time.
	 */
	Audio cut( Time start, Time end, Time startFade = 0, Time endFade = 0, flan_CANCEL_ARG ) const;

	/** This returns a peice of the original Audio that lied between start and end.
	 *	\param start Start frame for cut.
	 *	\param end End frame for cut.
	 *	\param startFade The cut audio start is faded for this amount of framea.
	 *	\param endFade The cut audio end is faded for this amount of frames.
	 */
	Audio cutFrames( Frame start, Frame end, Frame startFade = 0, Frame endFade = 0, flan_CANCEL_ARG ) const;

	/** This repitches the input.
	 *	\param factor The output of this represents a scaling factor of pitch as a function of time.
	 *	\param granularity This represents the frequency at which factor is sampled.
			Rapidly changing factors should use a small granularity.
	 *	\param qual This can take the value 0 for a default quality sinc resampling, 1 for quick and dirty (lerp and simple filter), 
	 *		or 2 for bad quality.
	 */
	Audio repitch( Func1x1 factor, Time granularity = .001f, uint32_t qual = 0, flan_CANCEL_ARG ) const;

	/** If each Func1x1 in ir is constant this is equivalent to normal convolution between the input Audio and ir.
	 *	\param ir An impulse response, this is treated like an Audio whos samples can change over time. 
	 */
	Audio convolve( const std::vector<Func1x1> & ir, flan_CANCEL_ARG ) const;

	/** This is normal convolution between Audio.
	 *	\param ir An impulse response.
	 */
	Audio convolve( const Audio & ir, flan_CANCEL_ARG ) const;

	/** This adds a fade to the ends of the input Audio.
	 *	\param start Length of the start fade.
	 *	\param end Length of the end fade.
	 *	\param interp The fading curve. Square root should be used for constant-power crossfading.
	 */
	Audio fade( Time start = 16.0f/48000.0f, Time end = 16.0f/48000.0f, 
		Interpolator interp = Interpolators::sqrt, flan_CANCEL_ARG ) const;

	/** This adds a fade to the ends of the input Audio.
	 *	\param start Length of the start fade.
	 *	\param end Length of the end fade.
	 *	\param interp The fading curve. Square root should be used for constant-power crossfading.
	 */
	Audio fadeFrames( Frame start = 16, Frame end = 16, 
		Interpolator interp = Interpolators::sqrt, flan_CANCEL_ARG ) const;

	/** A simple low pass filter. This method will likely be overhauled with the planned filter methods.
	 *	\param cutoff The frequency at which attenuation should begin as a function of time.
	 *	\param taps This in a sense describes the "accuracy" of the filter. 
	 *		Increasing it will produce a sharper cutoff and increase processing time.
	 */
	Audio lowPass( Func1x1 cutoff, uint32_t taps = 64, flan_CANCEL_ARG ) const;

	/** This repeats the input n times.
	 *	\param n The number of desired iterations
	 *	\param mod This allows user-defined processing of each iteration.
	 *	\param feedback When enabled, the previous iteration is sent to mod, rather than the original Audio.
	 */
	Audio iterate( uint32_t n, Mod mod = nullptr, bool feedback = false, flan_CANCEL_ARG ) const;

	/** Generates a volume-decaying iteration of the input. If a mod is supplied, it will 
	 *	be fed the previous iteration's ouput, rather than the original Audio.
	 *	\param length The time at which delay events stop generating.
	 *	\param delayTime The amount of time between delays.
	 *	\param decay Each delay event will have its gain scaled by the decay at the event time.
	 *	\param mod This allows user-defined processing of each delay.
	 */
	Audio delay( Time length, Func1x1 delayTime, Func1x1 decay = .5, 
		Audio::Mod mod = nullptr, flan_CANCEL_ARG ) const;

	/** Texture is the main Audio timing function. Events are generated following the
	 *	provided eventsPerSecond and scatter. Each event places a copy of the input Audio
	 *	at the event time. If a mod is provided, the Audio is first passed through the 
	 *	mod before placement.
	 *  \param length Events will generate until this time.
	 *  \param eventsPerSecond The mean number of events per second.
	 *  \param scatter The standard deviation in events. Due to being in events, a higher 
	 *		eventsPerSecond will cause scatter to have less of an effect.
	 *  \param mod Each event is fed into this along with the event time.
	 *  \param feedback This decides if mod will process the original Audio, or the previous mod output.
	 */
	Audio texture( Time length, Func1x1 eventsPerSecond, Func1x1 scatter = 0, Mod mod = nullptr, 
		bool feedback = false, flan_CANCEL_ARG ) const;

	/** This is a utility for calling Audio::texture with the most common mod functions. 
	 *  \param length Events will generate until this time.
	 *  \param eventsPerSecond The mean number of events per second.
	 *  \param scatter The standard deviation in events. Due to being in events, a higher 
	 *		eventsPerSecond will cause scatter to have less of an effect.
	 *  \param repitch Pitch scaling.
	 *  \param gain Volume scaling.
	 *	\param eventLength Event length.
	 *	\param pan Pan amount. Zero gives no panning.
	 */
	Audio textureSimple( Time length, Func1x1 eventsPerSecond, Func1x1 scatter, 
		Func1x1 repitch, Func1x1 gain, Func1x1 eventLength, Func1x1 pan = 0, flan_CANCEL_ARG ) const;

	/** This is a traditional granular synthesis tool.
	 *  \param length Grains will generate until this time.
	 *  \param grainsPerSecond The mean number of grains per second.
	 *  \param scatter The standard deviation in grains. Due to being in grains, a higher 
	 *		grainsPerSecond will cause scatter to have less of an effect.
	 *  \param selection The input time from which grains will be read.
	 *  \param grainLength The length of grains.
	 *  \param fade The start and end fade time for each grain. Fades use a sqrt curve.
	 *  \param mod Each grain is fed into this along with the event time.
	 */
	Audio grainSelect( Time length, Func1x1 grainsPerSecond, Func1x1 scatter, Func1x1 selection, 
		Func1x1 grainLength = 0.1, Func1x1 fade = 0.05, Mod mod = nullptr, flan_CANCEL_ARG ) const;


	/** This slices the input into a vector of chunks, each sliceLength seconds long.
	 *	The final chunk will be zero padded to meet the sliceLength.
	 *  \param sliceLength The length of each output.
	 *  \param fade The start and end fade time of each slice.
	 */
	std::vector<Audio> chop( Time sliceLength, Time fade = 0.05, flan_CANCEL_ARG ) const;

	/** This performs chop, randomizes the order of the chopped pieces, and joins, crossfading.
	 *  \param sliceLength The length of each output.
	 *  \param fade The start and end fade time of each slice.
	 */
	Audio rearrange( Time sliceLength, Time fade = 0.05, flan_CANCEL_ARG ) const;

	//========================================================
	// Multi-In Procs
	//========================================================

	/** The main Audio mixing method. The ouput format will be that of the first input.
	 *	\param ins The Audio to mix.
	 *	\param balances The volume scaling to apply to each input Audio. 
	 *		The default input will balance all inputs to 1/n, where n is the number of inputs.
	 *	\param startTimes This gives the start time for each input Audio. The default will start all inputs at time 0.
	 */
	static Audio mix( const std::vector<Audio> & ins, 
		std::vector<Func1x1> balances = std::vector<Func1x1>(),
		std::vector<Time> startTimes = std::vector<Time>(),
		flan_CANCEL_ARG );

	/** This joins all input Audio tip to tail in the order they were passed.
	 *	\param ins The Audio to mix.
	 */
	static Audio join( const std::vector<Audio> & ins, Time fade = 0, flan_CANCEL_ARG );

	/** At each point in time, selection decides which of the input Audio streams is playing.
	 *	Non-integer selections will mix appropriately scaled copies of the surrounding integer inputs.
	 *	\param ins The Audio to mix.
	 *	\param selection Which audio to play.
	 *	\param startTimes Time offsets for each input.
	 */
	static Audio select( const std::vector<Audio> & ins, Func1x1 selection,
		std::vector<Time> startTimes = std::vector<Time>(), flan_CANCEL_ARG );
};

} // End namespace flan
