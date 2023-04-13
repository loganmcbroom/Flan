#pragma once

#include <functional>
#include <complex>

#include "flan/Audio/AudioBuffer.h"
#include "flan/Audio/AudioMod.h"
#include "flan/Utility/Interpolator.h"
#include "flan/Utility/Interval.h"
#include "flan/Function.h"

namespace flan {

class Spectrum;
class PV;
class SPV;
class SQPV;
class Graph;
class FIRFilter;

//! Audio Algorithms
/** This is a wrapper of flan::AudioBuffer which contains all relevant algorithms and utilities available
 *	in flan for manipulating Audio buffers.
 *
 *	All Audio methods are const and have no side effects.
 */
class Audio : public AudioBuffer
{
public:

	Audio copy() const;
	
	/** Constructs 0 size Audio with default AudioBuffer::Format. 
	 */
	Audio();

	/** Constructs an audio from a temporary buffer.
	 */
	Audio( std::vector<float> && buffer, Channel numChannels, SampleRate );

	/** Constructs an Audio with the Format other but with an initialized buffer.
	 *	\param other Format for the constructed PV
	 */
	Audio( const AudioBuffer::Format & other );

	/** Constructs an Audio from the file at filename.
	 *	Supports any file extensions supported by libsndfile. Wave is suggested, and mp3 is unsupported.
	 *	\param filename Filename to load
	 */
	Audio( const std::string & filename );

	/** Audio is a stateless wrapper, so it can be constructed from a buffer
	 */
	Audio( AudioBuffer && );

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


	Audio resample( double newSampleRate ) const;

	/** Converts the Audio to a waveform bmp.
	 *  \param I The time interval to graph. Passing the defalt of (0,-1) will graph the entire Audio.
	 *  \param width The bmp width.
	 *  \param height The bmp height.
	 *	\param timelineScale Size of tick marks. Zero for none.
	 */
	Graph convertToGraph( Interval I = Interval( 0, -1 ), Pixel width = -1, Pixel height = -1, float timelineScale = 0 ) const;

	/** Converts the Audio to a waveform bmp and saves it at filename.
	 *	\param filename A filename at which to save the generated bmp.
	 *  \param I The interval of time to graph.
	 *  \param width The bmp width.
	 *  \param height The bmp height.
	 */
	void saveToBMP( const std::string & filename, Interval I = Interval( 0, -1 ), Pixel width = -1, Pixel height = -1 ) const;

	/** Apply a short-time Forier transform to the Audio, and phase vocode the output. See phase_vocoder for details on phase vocoding.
	 *
	 *	\param windowSize This is the number of Frames copied into each fft input. 
	 *		Defined behaviour is only guaranteed for power-of-two inputs, but it may work for other input sizes.
	 *		The number of bins in the output PV is given by frameSize / 2 + 1.
	 *	\param hop The number of Audio frames per PV frame.
	 *  \param dftSize The dft size. Unlike with the other time-frequency types that Audio can be converted to, this transform can use 
	 *		a dft size larger than the window.
	 */
	PV convertToPV( Frame windowSize = 2048, Frame hop = 128, Frame dftSize = 4096, flan_CANCEL_ARG ) const;

	/** For stereo inputs this is identical to Audio::convertToPV, but converts the audio to mid-side first. 
	 *	This is often better sounding than convertToPV because it removes the channel incoherece and phasing issues 
	 *	that multichannel PV processing often creates.
	 *	Non-stereo inputs will produce a null output.
	 *	See PV::convertToLeftRightAudio for the inverse transform.
	 */
	PV convertToMidSidePV( Frame windowSize = 2048, Frame hop = 128, Frame dftSize = 4096, flan_CANCEL_ARG ) const;

	/** Apply a sliding DFT to the Audio, and phase vocode the output. See phase_vocoder for details on phase vocoding. Be aware that this process
	 *	can return a very large output.
	 *
	 *  \param dftSize The dft size. This determines the number of frequency bins in the output.
	 */
	SPV convertToSPV( Frame dftSize = 1024 ) const;
	SPV convertToMidSideSPV( Frame dftSize = 1024 ) const;

	/** Apply a sliding constant-q transform to the Audio, and phase vocode the output. See phase_vocoder for details on phase vocoding. Be aware 
	 *	that this process can return a very large output.
     *
	 *  \param bandwidth The frequency range covered by the transform.
	 *	\param binsPerOctave The number of frequency bins per octave of bandwidth.
	 */
	SQPV convertToSQPV( std::pair<Frequency, Frequency> bandwidth = std::make_pair( 20, 22000 ), Bin binsPerOctave = 24 ) const;
	SQPV convertToMidSideSQPV( std::pair<Frequency, Frequency> bandwidth = std::make_pair( 20, 22000 ), Bin binsPerOctave = 24 ) const;

	/** This extracts a function approximating volume over time. The extracted function 
	 *	is generated by linear approximations. Only the first channel is processed.
	 *
	 *  \param granularity This parameter decides how often (in time) the current average gain is sampled. 
	 *		A granularity of 1 will return a function exactly equal to the input Audio.
	 *		Smaller granularities will generate more costly functions.
	 */
	Func1x1 convertToFunction( Time granularity = .001f ) const;

	/** AudioBuffer normally stores "Left" and "Right" data in channel 0 and 1 respectively.
	 *	The result of this transform stores "Mid" and "Side" data in those channels instead.
	 *	This allows some simple effects such as widening, but also allows mid and side information
	 *	to be processed independantly.
	 */
	Audio convertToMidSide() const;

	/** This converts back to Left/Right from Mid/Side. See Audio::convertToMidSide. 
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



	//============================================================================================================================================================
	// Information
	//============================================================================================================================================================

	/** Find the total energy in each channel of the input. Energy is the sum of the square of each sample. 
	 */
	std::vector<float> getTotalEnergy() const;

	/** Find the difference in energy between two Audios. This is useful for unit testing algorithms.
	 */
	std::vector<float> getEnergyDifference( const Audio & other ) const;

	/** Try to find the wavelength of the input over time. Selects to minimize differences per repitition.
	 *	\param channel The channel to process.
	 *	\param start Frame to analyze.
	 *	\param absoluteCutoff Minimum the local wave difference minimum must achieve to be considered a potential wavelength.
	 *	\param windowSize Input window size.
	 *	\param fft External fft handler. IMPORTANT: if using this, the fft size needs to be found with autocorrelationFFTSize( windowSize ) in DSPUtility.
	 */
	float getLocalWavelength( Channel channel, Frame start, Frame windowSize = 2048 ) const;

	/** Try to find the wavelength of the input over time. This is estimated using parabolic interpolation, so it may not be an integer.
	 *	\param channel The channel to process.
	 *	\param start Start frame.
	 *	\param end End frame.
	 *	\param absoluteCutoff Minimum the local wave difference minimum must achieve to be considered a potential wavelength.
	 *	\param windowSize The amount of input data used in each analysis.
	 *	\param hop The hop size in frames per anlysis.
	 */
	std::vector<float> getLocalWavelengths( Channel channel, Frame start = 0, Frame end = -1, Frame windowSize = 2048, 
		Frame hop = 128, flan_CANCEL_ARG ) const;

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
	 */
	float getAverageWavelength( Channel, float minActiveRatio = 0, float maxLengthSigma = -1, Frame start = 0, Frame end = -1, 
		Frame windowSize = 2048, Frame hop = 128, flan_CANCEL_ARG ) const; // Averages getLocalWavelengths
	

	/** Try to find the wavelength of the input over time. Selects highest to minimize octave error.
	 *	\param channel The channel to process.
	 *	\param start Frame to analyze.
	 *	\param windowSize Input window size.
	 *	\param fft External fft handler. IMPORTANT: if using this, the fft size needs to be found with autocorrelationFFTSize( windowSize ) in DSPUtility.
	 */
	Frequency getLocalFrequency( Channel channel, Frame start = 0, Frame windowSize = 2048, flan_CANCEL_ARG ) const;

	/** Utility for calling getLocalFrequency at a set frame hop interval from start to end.
	 */
	std::vector<Frequency> getLocalFrequencies( Channel channel, Frame start = 0, Frame end = -1, Frame windowSize = 2048, 
		Frame hop = 128, flan_CANCEL_ARG ) const;

	//============================================================================================================================================================
	// Processing
	//============================================================================================================================================================

	/** This phase inverts the input. Every output sample is assigned the negative of the input sample at that time.
	 */
	Audio invertPhase() const;

	/** Volume scaling, Output( t ) = input( t ) * volumeLevel( t ).
	 *	\param gain Describes how volume should be scaled as a function of time.
	 */
	Audio modifyVolume( const Func1x1 & gain ) const;
	void modifyVolumeInPlace( const Func1x1 & gain );

	/** Volume setting. This normalizes (divides every sample by the largest sample value), and then scales by level.
	 *	\param level Describes how volume should be scaled as a function of time.
	 */
	Audio setVolume( const Func1x1 & level ) const;

	/** Time shifting.
	 *	\param shift How much to shift the input.
	 */
	Audio shift( Time shift ) const;

	/** This applies the shaper as a function to each sample in the input.
	 *	\param shaper Each sample in the input is passed through this. 
	 *		Samples will be values on [-1,1] under normal circumstances.
	 *		For example, the function y = x would have no effect as a shaper.
	 */
	Audio waveshape( const Func1x1 & shaper ) const;

	/// This can't be implemented until I have a setup with more than two speakers to test it on
	// /** This spatially repositions the input. This is a general purpose panning algorithm for handling any number of speakers. 
	//  *	The speakers are assumed to be in a plane. Units are given in meters for haas effect purposes.
	//  *	\param panPosition The plane position of the input over time. 
	//  *  \param listenerPosition The plane position of the listener over time.
	//  *	\param speakerPositions The plane positions of each speaker.
	//  *  \param hassEffect Decides if the haas effect should be applied. 
	//  */
	// Audio pan( Function<Time,vec2> panPosition, Function<Time, vec2> listenerPosition, std::vector<vec2> speakerPositions, bool haasEffect = false ) const;

	/** This spatially repositions the input for mono and stereo inputs. Stereo speakers are assumed form an equilateral trangle with side
	 *	length 1 meter with the listener. The input position is assumed to sit halfway between the speakers.
	 *  If the number of channels is neither 1 nor 2 a null Audio is returned.
	 *	\param panPosition The input spatial position from -1 (left) to 1 (right) over time.
	 */
	Audio pan( const Func1x1 & panPosition ) const;
	void panInPlace( const Func1x1 & panPosition );

	/** This redistributes energy between the mid and side signals
	 *	\param widenAmount This should return a value on [-1,1], representing movement from mid to side.
	 */
	Audio widen( const Func1x1 & widenAmount ) const;

	/** This reverses the Audio.
	 */
	Audio reverse() const;

	/** This returns a peice of the original Audio that lied between startTime and endTime.
	 *	\param start Start time for cut.
	 *	\param end End time for cut.
	 *	\param startFade The cut audio start is faded for this amount of time.
	 *	\param endFade The cut audio end is faded for this amount of time.
	 */
	Audio cut( Time start, Time end, Time startFade = 0, Time endFade = 0 ) const;

	/** This returns a peice of the original Audio that lied between start and end.
	 *	\param start Start frame for cut.
	 *	\param end End frame for cut.
	 *	\param startFade The cut audio start is faded for this amount of framea.
	 *	\param endFade The cut audio end is faded for this amount of frames.
	 */
	Audio cutFrames( Frame start, Frame end, Frame startFade = 0, Frame endFade = 0 ) const;

	enum class WDLResampleType
		{
		Sinc = 0,
		Linear = 1,
		Uninterpolated = 2,
		};
	/** This repitches the input.
	 *	\param factor The output of this represents a scaling factor of pitch as a function of time.
	 *	\param granularity This represents the frequency at which factor is sampled.
			Rapidly changing factors should use a small granularity.
	 *	\param qual This can take the value 0 for a default quality sinc resampling, 1 for quick and dirty (lerp and simple filter), 
	 *		or 2 for bad quality.
	 */
	Audio repitch( const Func1x1 & factor, Time granularity = .001f, WDLResampleType qual = WDLResampleType::Sinc ) const;

	/** If each Func1x1 in ir is constant this is equivalent to normal convolution between the input Audio and ir.
	 *	\param ir An impulse response, this is treated like an Audio whos samples can change over time. 
	 */
	Audio convolve( const Function<Time, std::vector<float>> & ir ) const;

	/** This is normal convolution between Audio. If the number of channels are mismatched, the ir channels will be used cyclically.
	 *	\param ir An impulse response.
	 */
	Audio convolve( const Audio & ir ) const;

	/** This adds a fade to the ends of the input Audio.
	 *	\param start Length of the start fade.
	 *	\param end Length of the end fade.
	 *	\param interp The fading curve. Square root should be used for constant-power crossfading.
	 */
	Audio fade( Time start = 16.0f/48000.0f, Time end = 16.0f/48000.0f, Interpolator interp = Interpolators::sqrt ) const;
	void fadeInPlace( Time start = 16.0f/48000.0f, Time end = 16.0f/48000.0f, Interpolator interp = Interpolators::sqrt );

	/** This adds a fade to the ends of the input Audio.
	 *	\param start Length of the start fade.
	 *	\param end Length of the end fade.
	 *	\param interp The fading curve. Square root should be used for constant-power crossfading.
	 */
	Audio fadeFrames( Frame start = 16, Frame end = 16, Interpolator interp = Interpolators::sqrt ) const;
	void fadeFramesInPlace( Frame start = 16, Frame end = 16, Interpolator interp = Interpolators::sqrt );

	// /** A simple low pass filter. This method will likely be overhauled with the planned filter methods.
	//  *	\param cutoff The frequency at which attenuation should begin as a function of time.
	//  *	\param taps This in a sense describes the "accuracy" of the filter. 
	//  *		Increasing it will produce a sharper cutoff and increase processing time.
	//  */
	// Audio lowPass( const Func1x1 & cutoff, uint32_t taps = 64, flan_CANCEL_ARG ) const;

	/** This repeats the input n times.
	 *	\param n The number of desired iterations
	 *	\param mod This allows user-defined processing of each iteration.
	 *	\param feedback When enabled, the previous iteration is sent to mod, rather than the original Audio.
	 */
	Audio iterate( uint32_t n, const AudioMod & mod = AudioMod(), bool feedback = false ) const;

	/** Generates a volume-decaying iteration of the input. If a mod is supplied, it will 
	 *	be fed the previous iteration's ouput, rather than the original Audio.
	 *	\param length The time at which delay events stop generating.
	 *	\param delayTime The amount of time between delays.
	 *	\param decay Each delay event will have its gain scaled by the decay at the event time.
	 *	\param mod This allows user-defined processing of each delay.
	 */
	Audio delay( Time length, const Func1x1 & delayTime, const Func1x1 & decay = .5, const AudioMod & mod = AudioMod() ) const;

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
	Audio texture( Time length, const Func1x1 & eventsPerSecond, const Func1x1 & scatter = 0, const AudioMod & mod = AudioMod(), bool feedback = false ) const;

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
	Audio textureSimple( Time length, const Func1x1 & eventsPerSecond, const Func1x1 & scatter, 
		const Func1x1 & repitch, const Func1x1 & gain, const Func1x1 & eventLength, const Func1x1 & pan = 0 ) const;

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
	Audio grainSelect( Time length, const Func1x1 & grainsPerSecond, const Func1x1 & scatter, const Func1x1 & selection, 
		const Func1x1 & grainLength = 0.1, const Func1x1 & fade = 0.05, const AudioMod & mod = AudioMod() ) const;


	/** This slices the input into a vector of chunks, each sliceLength seconds long.
	 *	The final chunk will be zero padded to meet the sliceLength.
	 *  \param sliceLength The length of each output.
	 *  \param fade The start and end fade time of each slice.
	 */
	std::vector<Audio> chop( Time sliceLength, Time fade = 0.05 ) const;

	/** This performs chop, randomizes the order of the chopped pieces, and joins, crossfading.
	 *  \param sliceLength The length of each output.
	 *  \param fade The start and end fade time of each slice.
	 */
	Audio rearrange( Time sliceLength, Time fade = 0.05 ) const;

	//========================================================
	// Static Procs
	//========================================================

	/** The main Audio mixing method. The ouput format will be that of the first input.
	 *	\param ins The Audio to mix.
	 *	\param startTimes This gives the start time for each input Audio. The default will start all inputs at time 0.
	 *	\param balances The volume scaling to apply to each input Audio. The default will apply no scaling.
	 *		These are relative to global time 0 rather than the start time of each input.
	 */
	static Audio mix( std::vector<const Audio *> ins, 
					  std::vector<Time> startTimes = {}, 
					  const std::vector<const Func1x1 *> & gains = {} );

	static Audio mix( const std::vector<Audio> & ins, 
					  const std::vector<Time> & startTimes = {}, 
					  const std::vector<Func1x1> & gains = {} );

	/** This joins all input Audio tip to tail in the order they were passed.
	 *	\param ins The Audio to join.
	 *	\param offset The amount of time away from the natural join position each input should be. For example, using a negative offset
	 *		along with fading the inputs can be used to crossfade.
	 */
	static Audio join( const std::vector<const Audio *> & ins, Time offset = 0 );
	static Audio join( const std::vector<Audio> & ins, Time offset = 0 );

	/** At each point in time, selection decides which of the input Audio streams is playing.
	 *	Non-integer selections will mix appropriately scaled copies of the surrounding integer inputs.
	 *	\param ins The Audio to mix.
	 *	\param selection Which audio to play.
	 *	\param startTimes Time offsets for each input.
	 */
	static Audio select( const std::vector<const Audio *> & ins, const Func1x1 & selection, const std::vector<Time> & startTimes = {} );
	 static Audio select( const std::vector<Audio> & ins, const Func1x1 & selection, std::vector<Time> startTimes = {} );

	/** Generate an Audio from a Func1x1. This can be accomplished in a more general setting using the function ctor of flan::Wavetable.
	*  \param wave This is evaluated from 0 to 1. This portion is repeated to create the Audio.
	*  \param length The length of the output.
	*  \param freq The frequency of the output.
	*  \param sampleRate The sample rate of the output.
	*  \param oversample Synthesized audio must be generated at a sample rate higher than the desired sample rate to avoid aliasing.
	*      This describes how many samples should be used in the synthesis per sample in the output.
	*/
	static Audio synthesize( const Func1x1 & wave, Time length, const Func1x1 & freq, size_t samplerate = 48000, size_t oversample = 16 );

	};

} // End namespace flan
