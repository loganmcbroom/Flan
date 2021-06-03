#pragma once

#include <functional>
#include <map>

#include "PVOCBuffer.h"
#include "xcdp/Function.h"
#include "xcdp/Utility/Interpolator.h"
#include "xcdp/Utility/Rect.h"

struct fftwf_plan_s;
namespace std { template<class T> class complex; };

namespace xcdp {

class Audio;
class FFTHelper;
class Graph;

//! Phase Vocoder (PVOC) Data Algorithms
/** This is a stateless wrapper of xcdp::PVOCBuffer which contains all relevant algorithms and utilities available
 *	in xcdp for manipulating PVOC buffers.
 *
 *	All methods are const and have no side effects outside of PVOC::graph, which saves a bmp file.
 *	See PVOCBuffer for a description of Phase Vocoder (PVOC) Data.
 */
class PVOC : public PVOCBuffer
{
public:

	using Vec = const std::vector<const PVOC &> &;

	/** Constructs 0 size PVOC with default PVOCBuffer::Format. */
	PVOC() : PVOCBuffer( PVOCBuffer::Format() ) {}

	/** Constructs a PVOC with the Format other but with an initialized buffer.
	 *	\param other Format for the constructed PVOC
	 */
	PVOC( const PVOCBuffer::Format & other ) : PVOCBuffer( other ) {}

	/** Constructs a PVOC from the .ana file at filename.
	 *	\param filename A .pvoc filename to load. 
			This file type is xcdp specific and is generated by PVOCBuffer::save.
	 */
	PVOC( const std::string & filename ) : PVOCBuffer( filename ) {}

	/** PVOC is a stateless wrapper, so it can be constructed from a buffer
	 */
	PVOC( const PVOCBuffer & other ) : PVOCBuffer( other ) {}



	//============================================================================================================================================================
	// Conversions
	//============================================================================================================================================================

	/** Transforms this PVOC to an Audio using the overlaps and samplerate given in format. 
	 *  \param fftwRealBuffer When fftw is used this allows allocating the real fftw buffer used in the transform
	 *		on the main thread, as fftw buffer allocations aren't thread safe.
	 *  \param fftwOutBuffer When fftw is used this allows allocating the complex fftw buffer used in the transform
	 *		on the main thread, as fftw buffer allocations aren't thread safe.
	 *	 \param plan When fftw is used this allows allocating the fftw plan used in the transform
	 *		on the main thread, as fftw planning isn't thread safe.
	 */
	Audio convertToAudio( std::shared_ptr<FFTHelper> fft = nullptr, XCDP_CANCEL_ARG ) const;

	/** See PVOC::convertToAudio.
	 */ 
	Audio convertToAudio_cpu( std::shared_ptr<FFTHelper> fft = nullptr, XCDP_CANCEL_ARG ) const;

	/** Converts the PVOC to a spectrograph bmp.
	 *  \param domain The time/frequency rectangle to graph. Negative 1 for time or frequency end will use the maximum.
	 *  \param width The bmp width.
	 *  \param height The bmp height.
	 */
	Graph convertToGraph( Rect domain = { 0, 0, -1, -1 }, Pixel width = -1, Pixel height = -1, float timelineScale = 0, XCDP_CANCEL_ARG ) const;

	/** Creates and saves a bmp spectrograph of the PVOC, then returns this.
	 *	\param filename A filename at which to save the generated bmp.
	 *  \param domain The time/frequency rectangle to graph. Negative 1 for time or frequency end will use the maximum.
	 *  \param width The bmp width.
	 *  \param height The bmp height.
	 */
	const PVOC & saveToBMP( const std::string & filename, Rect domain = { 0, 0, -1, -1 }, Pixel width = -1, Pixel height = -1, XCDP_CANCEL_ARG ) const;



	//============================================================================================================================================================
	// Contours
	//============================================================================================================================================================

	/** Container for salience (percieved volume) extracted by PVOC::getSalience.
	 */
	struct Salience
		{
		float & get( Frame f, Bin b ) { return buffer[ f * numBins + b ]; }
		Frame numFrames;
		Bin numBins;
		std::vector<float> buffer;
		};

	/** This attempts to find the percieved volume at every tenth of a note over time. Meant to be used by PVOC::getContours.
	 *    Note this doesn't use a percieved volume filter like it maybe should. I didn't want to deal with it.
	 *	\param channel The channel to analyze.
	 *	\param minFreq Frequency of the lowest pitch bin.
	 *	\param maxFreq Frequency of the highest pitch bin.
	 */
	Salience getSalience( Channel channel, Frequency minFreq = 55, Frequency maxFreq = 1760, XCDP_CANCEL_ARG ) const;

	/** See PVOC::getSalience.
	 */
	Salience getSalience_gpu( Channel channel, Frequency minFreq = 55, Frequency maxFreq = 1760, XCDP_CANCEL_ARG ) const;

	/** Container for a contour (note) extracted by PVOC::getContours, with some commonly needed details.
	 */
	struct Contour
		{
		float pitchMean;
		float pitchSD;
		float salienceMean;
		float salienceSD;

		Frame startFrame;
		std::vector<vec2> bins; // Each vec2 contains { pitchBin, mag (?) }
		};

	/** This is a piece of the melodia melody finding algorithm, which attempts to extract note information from PVOC data.
	 *  \param channel The channel to analyze.
	 *  \param minFreq
	 *  \param maxFreq
	 *  \param filterShort Contours with lengths less than this will be removed.
	 *  \param filterQuiet Contours with salienceMean less than the largest salienceMean divided by this are removed
	 */
	std::vector<Contour> getContours( Channel channel, Frequency minFreq = 55, Frequency maxFreq = 1760, Frame filterShort = 30, float filterQuiet = 20, XCDP_CANCEL_ARG ) const;

	/** Functuon type exclusive to PVOC::prism. */
	using PrismFunc = std::function<MF ( Time, int harmonic, Frequency contourFreq, const std::vector<Magnitude> & harmonicMagnitudes )>;

	/** Prism, in theory, allows complete control over the frequency and magnitude of every harmonic of every note in the input.
	 *  \param harmonicFunction This takes the global or local time (see perNote), the harmonic index (starting at 1), the base frequency and the magnitudes of all harmonics
	 *  \param perNote This decides if the time passed to harmonicFunction should be the time elapsed since the start of the PVOC (false), or the start of the contour (true).
	 */
	PVOC prism( PrismFunc harmonicFunction, bool perNote = true, XCDP_CANCEL_ARG ) const;

	//============================================================================================================================================================
	// Utility
	//============================================================================================================================================================

	/** Returns a single frame from the input. 
	 *	The returned frame in a linear interpolation of the frame surrounding the selected time.
	 *	\param time The time at which the frame should be taken. 
	 */
	PVOC getFrame( Time time ) const;

	/** This computes a weighted approximation of the 4 surrounding bins using the provided interpolator */
	MF getBinInterpolated( Channel channel, float frame, float bin, Interpolator interp = Interpolators::linear ) const;
	MF getBinInterpolated( Channel channel, float frame, Bin bin, Interpolator interp = Interpolators::linear ) const;
	MF getBinInterpolated( Channel channel, Frame frame, float bin, Interpolator interp = Interpolators::linear ) const;



	//============================================================================================================================================================
	// Selection
	//============================================================================================================================================================

	/** Every time/frequency point in the output PVOC will read an arbitrary point in the input determined by selector.
	 *	Selected non-integer coordinates will interpolate surrounding integer coordinates using interp.
	 *	\param length The length of the output PVOC
	 *	\param selector This function takes each time/frequency point in the output and returns 
	 *		the time/frequency that should be read from the input into that output
	 *	\param interp This determines how input points surrounding selected points are interpolated
	 */
	PVOC select( Time length, Func2x2 selector, Interpolator interp = Interpolators::linear, XCDP_CANCEL_ARG ) const;

	/** This is a classic "time-freeze" effect. At supplied times playback is "frozen", repeating the current frame 
	 *	(or a linear interpolation of surrounding frames), for a specified amount of time. After this time has elapsed,
	 *	playback resumes until another freeze time is encountered.
	 *	\param timing Each element of timing describes a point of time to freeze, and for 
	 *		how long it should be frozen, in that order
	 */
	PVOC freeze( const std::vector<std::array<Time,2>> & timing, XCDP_CANCEL_ARG ) const;



	//============================================================================================================================================================
	// Resampling
	//============================================================================================================================================================

	/** Each input bin is mapped to an arbitrary output point by mod. 
	 *	Every quad in the input is thus mapped to a quad in the output. 
	 *	The output points within that quad are given an average of any data already within the bin 
	 *	and the four input bins that made up the quad being mapped.
	 *	Quad averages are decided using the interpolator argument and the algorithm described here: 
	 *	https://www.particleincell.com/2012/quad-interpolation/
	 *	\param mod Takes and returns time/frequency pairs as either xcdp::vec2 or std::complex<float>
	 *	\param interp Interpolator used in quad mapping
	 */
	PVOC modify( Func2x2 mod, Interpolator interp = Interpolators::linear, XCDP_CANCEL_ARG ) const;

	/** This is functionally equivalent to using PVOC::modify and only outputting the input time
	 *	\param mod Takes time/frequency pairs and returns frequency
	 *	\param interp Interpolator used in frequency mapping
	 */
	PVOC modifyFrequency( Func2x1 mod, Interpolator = Interpolators::linear, XCDP_CANCEL_ARG ) const;

	/** This is functionally equivalent to using PVOC::modify and only outputting the input frequency
	 *	\param mod Takes time/frequency pairs and returns time
	 *	\param interp Interpolator used in time mapping
	 */
	PVOC modifyTime( Func2x1 mod, Interpolator = Interpolators::linear, XCDP_CANCEL_ARG ) const;

	/** This is nearly identical to PVOC::modify.
	 *	Due to the computations being on the gpu in parallel it is much faster.
	 *	The cost of this speed is that the output buffer is write only, so any output bin mapped to more than once
	 *  will take one of those written values only, with no determination as to which.
	 *	This method will return the input if xcdp was built without OpenCL.
	 *	\param mod Takes and returns time/frequency pairs as either xcdp::vec2 or std::complex<float>
	 *	\param interp Disabled for gpu methods.
	 */
	PVOC modify_cl( Func2x2 mod, Interpolator interp = Interpolators::linear, XCDP_CANCEL_ARG ) const;

	/** This is functionally equivalent to using PVOC::modify_cl and only outputting the input time
	 *	See PVOC::modify_cl for the effects of non-parallel computation
	 *	This method will return the input if xcdp was built without OpenCL.
	 *	\param mod Takes time/frequency pairs and returns frequency
	 *	\param interp Disabled for gpu methods.
	 */
	PVOC modifyFrequency_cl( Func2x1 mod, Interpolator interp = Interpolators::linear, XCDP_CANCEL_ARG ) const;

	/** This is functionally equivalent to using PVOC::modify_cl and only outputting the input frequency
	 *	See PVOC::modify_cl for the effects of non-parallel computation
	 *	This method will return the input if xcdp was built without OpenCL.
	 *	\param mod Takes time/frequency pairs and returns time
	 *	\param interp Interpolator used in time mapping, this is sampled 512 times for use by OpenCL
	 */
	PVOC modifyTime_cl( Func2x1 mod, Interpolator interp = Interpolators::linear, XCDP_CANCEL_ARG ) const;

	/** This is functionally equivalent to using PVOC::modifyFrequency_cl with the mod output multiplied by input frequency
	 *	\param factor Takes time/frequency pairs and returns frequency multiplier
	 *	\param interp Disabled for gpu methods.
	 */
	PVOC repitch( Func2x1 factor, Interpolator = Interpolators::linear, XCDP_CANCEL_ARG ) const;

	/** This is functionally equivalent to using PVOC::modifyTime_cl with the mod output multiplied by input time
	 *	\param factor Takes time/frequency pairs and returns time multiplier
	 *	\param interp Interpolator used in frequency mapping
	 */
	PVOC stretch( Func2x1 factor, Interpolator = Interpolators::linear, XCDP_CANCEL_ARG ) const;

	/** This is close to PVOC::stretch, but can only expand the input by integer quantities at any given time.
	 *	It is also restricted to only choosing the local expansion amount as a function of time.
	 *	Spline interpolation is used to fill the expanded space.
	 *	\param interp This returns the local expansion amount as a function of time. 
	 *		Output will be rounded to the nearest integer and clamped to positive integers.
	 */
	PVOC stretch_spline( Func1x1 interp, XCDP_CANCEL_ARG ) const;

	/** This is functionally equivalent to stretching the input down by factor, and then up.
	 *	Resolution is lost, but the lost resolution is filled via interpolation
	 *	\param factor Local information loss ratio. A factor of two will replace every other frame with interpolated data.
	 *	\param interp Interpolator deciding how lost frames are restored from their surroundings
	 */
	PVOC desample( Func2x1 factor, Interpolator interp = Interpolators::linear, XCDP_CANCEL_ARG ) const;

	/** Warning: this process can produce very loud outputs with unscrupulouss parameters.
	 *	The input is left alone until startTime. Between startTime and endTime the output is an interpolation
	 *	of the frames at startTime and endTime. After endTime, the interpolation continues into the realm of extrapolation
	 *	for a duration of extrapDuration.
	 *	\param startTime Time at which interpolation starts.
	 *	\param endTime Time at which interpolation ends and extrapolation begins
	 *	\param extrapDuration Duration of the extrapolation
	 *	\param interp Interpolator used throughout
	 */
	PVOC timeExtrapolate( Time startTime, Time endTime, Time extrapDuration, 
		Interpolator interp = Interpolators::linear, XCDP_CANCEL_ARG ) const;



	//============================================================================================================================================================
	// Extras
	//============================================================================================================================================================

	/** For all frequencies, every octave above it is set to a copy of the base, scaled by seriesScale
	 *	\param seriesScale The scaling function. The inputs to this are time and harmonic index starting at 0.
	 */
	PVOC addOctaves( Func2x1 seriesScale, XCDP_CANCEL_ARG ) const;

	/** For all frequencies, every harmonic above it is set to a copy of the base, scaled by seriesScale
	 *	\param seriesScale The scaling function. The inputs to this are time and harmonic index starting at 0.
	 */
	PVOC addHarmonics( Func2x1 seriesScale, XCDP_CANCEL_ARG ) const;

	/** Replaces the Amplitudes (Magnitudes) of bins in the input with those in ampSource
	 *	\param ampSource Source PVOC from which to draw amplitudes
	 *	\param amount An amount of 1 fully replaces the amplitudes, 0 does nothing, and amounts between
	 *		give a linear interpolation of the two.
	 *	\returns A PVOC with dimensions equal to the shorter of this and ampSource
	 */
	PVOC replaceAmplitudes( const PVOC & ampSource, Func2x1 amount = 1, XCDP_CANCEL_ARG ) const;

	/** Subtracts the Amplitudes (Magnitudes) of bins in other with those in this
	 *	\param ampSource Source PVOC from which to draw amplitudes to subtract.
	 *	\param amount A scalar on the subtraction amount. This isn't clamped to [0,1].
	 *	\returns A PVOC with dimensions equal to the shorter of this and ampSource.
	 */
	PVOC subtractAmplitudes( const PVOC & other, Func2x1 amount = 1, XCDP_CANCEL_ARG ) const;

	/** Passes each MF through the shaper.
	 *	\param shaper A function taking MFs and returning MFs.
	 */
	PVOC shape( Func2x2 shaper, XCDP_CANCEL_ARG ) const;

	/** Modifies the input data in random ways using normal distributions.
	 *	\param magSigma Standard deviation for magnitude distribution.
	 *	\param frqSigma Standard deviation for frequency distribution.
	 */
	PVOC perturb( Func2x1 magSigma, Func2x1 frqSigma, XCDP_CANCEL_ARG ) const;

	/** At any given time, numPartials should return a number of bins, N, to retain.
	 *	The N loudest bins are copied to the output.
	 *	All other output bins are 0 filled.
	 *	\param numPartials Number of partials to retain as a function of time
	 */
	PVOC retainNLoudestPartials( Func1x1 numPartials, XCDP_CANCEL_ARG ) const;

	/** At any given time, numPartials should return a number of bins, N, to remove.
	 *	The N loudest bin positions are 0 filled in the output.
	 *	All other bins are copied to the output.
	 *	\param numPartials Number of partials to remove as a function of time
	 */
	PVOC removeNLoudestPartials( Func1x1 numPartials, XCDP_CANCEL_ARG ) const;

	/** Each bin has its frequency copied to subsequent output bins with decaying magnitude until a 
	 *	bin with magnitude greater than the current decayed magnitude is read from the input. 
	 *	The frequency of the output at that point then changes to that of the louder bin.
	 *	\param decay This returns the decay amount per second at every time/frequency input.
	 *		For example, a constant decay of .5 applied to an impulse will lose half it's magnitude every second.
	 *	\param length Because the decay is exponential, an output length is needed. 
	 */
	PVOC resonate( Func2x1 decay, Time length, XCDP_CANCEL_ARG ) const;



	//============================================================================================================================================================
	// Backend
	//============================================================================================================================================================


};

} // End namespace xcdp

