#pragma once

#include "flan/Audio/AudioBuffer.h"
#include "flan/Audio/AudioMod.h"
#include "flan/Utility/Interpolator.h"
#include "flan/Utility/Interval.h"
#include "flan/Utility/vec2.h"
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
 *	Nearly all Audio methods are const and have no side effects. Impure functions have an "_in_place" suffix.
 */
class Audio : public AudioBuffer
{
public:

	/** SFINAE type for variadic templates, see https://www.fluentcpp.com/2019/01/25/variadic-number-function-parameters-type/
	 */
	template<typename... Ts>
	using AllAudios = typename std::enable_if_t<std::conjunction_v< std::is_convertible<Ts, const Audio &>... >>;

	template<typename O>
	FunctionSample<O> sample_function_over_domain( const Function<Second, O> & f ) const	
		{
		return f.sample( 0, get_num_frames(), 1.0f / get_sample_rate() );
		}

	//============================================================================================================================================================
	// Constructors
	//============================================================================================================================================================
	
	Audio();

	/** Audio is a stateless wrapper, so it can be constructed from a buffer
	 */
	Audio( 
		AudioBuffer && 
		);

	Audio copy(
		) const;
	
	/** Constructs 0 size Audio with default AudioBuffer::Format. 
	 */
	static Audio create_null(
		);

	/** Constructs an audio from a temporary buffer.
	 */
	static Audio create_from_buffer( 
		std::vector<float> && buffer, 
		Channel num_channels, 
		FrameRate sample_rate 
		);

	/** Constructs an Audio with the Format other but with an initialized buffer.
	 *	\param other Format for the constructed PV
	 */
	static Audio create_from_format( 
		const AudioBuffer::Format & other 
		);

	/** Constructs an Audio from the file at filename.
	 *	Supports any file extensions supported by libsndfile. Wave is suggested, and mp3 is unsupported.
	 *	\param filename Filename to load
	 */
	static Audio load_from_file( 
		const std::string & filename 
		);
	static Audio load_from_file( 
		const std::string & filename,
		SndfileStrings &
		);

	static Audio create_empty_with_length( 
		Second length, 
		Channel num_channels = 1, 
		FrameRate sample_rate = 48000.0f 
		);

	static Audio create_empty_with_frames( 
		Frame num_frames, 
		Channel num_channels = 1, 
		FrameRate sample_rate = 48000.0f 
		);



	//============================================================================================================================================================
	// Conversions
	//============================================================================================================================================================

	Audio resample( 
		FrameRate new_sample_rate 
		) const; 

	/** Converts the Audio to a waveform bmp.
	 *  \param I The time interval to graph. Passing the defalt of (0,-1) will graph the entire Audio.
	 *  \param width The bmp width.
	 *  \param height The bmp height.
	 *	\param timeline_scale Size of tick marks. Zero for none.
	 */
	Graph convert_to_graph( 
		Interval I = Interval( 0, -1 ), 
		Pixel width = -1, 
		Pixel height = -1, 
		Graph::WaveformMode mode = Graph::WaveformMode::Symmetric,
		float timeline_scale = 0 
		) const;

	/** Converts the Audio to a waveform bmp and saves it at filename.
	 *	\param filename A filename at which to save the generated bmp.
	 *  \param I The interval of time to graph.
	 *  \param width The bmp width.
	 *  \param height The bmp height.
	 */
	void save_to_bmp( 
		const std::string & filename, 
		Interval I = Interval( 0, -1 ), 
		Pixel width = -1, 
		Pixel height = -1 
		) const;

	Graph convert_to_spectrum_graph(
		Pixel width = -1, 
		Pixel height = -1,
		Frame smoothing_frames = 128
		) const; 

	void save_spectrum_to_bmp( 
		const std::string & filename, 
		Pixel width = -1, 
		Pixel height = -1,
		Frame smoothing_frames = 128
		) const;

	/** Apply a short-time Forier transform to the Audio, and phase vocode the output. See phase_vocoder for details on phase vocoding.
	 *
	 *	\param window_size This is the number of Frames copied into each fft input. 
	 *		Defined behaviour is only guaranteed for power-of-two inputs, but it may work for other input sizes.
	 *		The number of bins in the output PV is given by frameSize / 2 + 1.
	 *	\param hop The number of Audio frames per PV frame.
	 *  \param dft_size The dft size. Unlike with the other time-frequency types that Audio can be converted to, this transform can use 
	 *		a dft size larger than the window.
	 */
	PV convert_to_PV( 
		Frame window_size = 2048, 
		Frame hop = 128, 
		Frame dft_size = 4096, 
		flan_CANCEL_ARG 
		) const;

	/** For stereo inputs this is identical to Audio::convert_to_PV, but converts the audio to mid-side first. 
	 *	This is often better sounding than convert_to_PV because it removes the channel incoherece and phasing issues 
	 *	that multichannel PV processing often creates.
	 *	Non-stereo inputs will produce a null output.
	 *	See PV::convert_to_lr_audio for the inverse transform.
	 */
	PV convert_to_ms_PV( 
		Frame window_size = 2048, 
		Frame hop = 128, 
		Frame dft_size = 4096, 
		flan_CANCEL_ARG 
		) const;

	/** Apply a sliding DFT to the Audio, and phase vocode the output. See phase_vocoder for details on phase vocoding. Be aware that this process
	 *	can return a very large output.
	 *
	 *  \param dft_size The dft size. This determines the number of frequency bins in the output.
	 */
	SPV convert_to_SPV( 
		Frame dft_size = 1024 
		) const;

	SPV convert_to_ms_SPV( 
		Frame dft_size = 1024 
		) const;

	/** Apply a sliding constant-q transform to the Audio, and phase vocode the output. See phase_vocoder for details on phase vocoding. Be aware 
	 *	that this process can return a very large output.
     *
	 *  \param bandwidth The frequency range covered by the transform.
	 *	\param bins_per_octave The number of frequency bins per octave of bandwidth.
	 */
	SQPV convert_to_SQPV( 
		std::pair<Frequency, Frequency> bandwidth = std::make_pair( 16, 24000 ), 
		Bin bins_per_octave = 24 
		) const;

	SQPV convert_to_ms_SQPV( 
		std::pair<Frequency, Frequency> bandwidth = std::make_pair( 16, 24000 ), 
		Bin bins_per_octave = 24 
		) const;

	/** AudioBuffer normally stores "Left" and "Right" data in channel 0 and 1 respectively.
	 *	The result of this transform stores "Mid" and "Side" data in those channels instead.
	 *	This allows some simple effects such as widening, but also allows mid and side information
	 *	to be processed independantly.
	 */
	Audio convert_to_mid_side(
		) const;

	/** This converts back to Left/Right from Mid/Side. See Audio::convert_to_mid_side. 
	 */
	Audio convert_to_left_right(
		) const;

	/** This converts the input to a two channel Audio. The conversion is channel-count dependant.
	 *	The currently accepted channel counts are: 1 and 2.
	 */
	Audio convert_to_stereo(
		) const;

	/** This converts an Audio with any number of channels to a single channel Audio.
	 *	Each output sample is the sum of the channels at that time over the square root of the number of channels.
	 */
	Audio convert_to_mono(
		) const;

	Function<Second, Amplitude> convert_to_function(
		) const;



	//============================================================================================================================================================
	// Channels
	//============================================================================================================================================================

	std::vector<Audio> split_channels(
		) const;

	static Audio combine_channels( 
		const std::vector<const Audio *> & channels 
		);

	static Audio combine_channels( 
		const std::vector<Audio> & channels 
		);

	/** Creates a new Audio with all the channels of the inputs.
	 */
	template<typename ... Ts, typename = AllAudios<Ts...> >
	static Audio combine_channels( 
		const Ts & ... ins 
		)
		{
		std::vector<const Audio *> p_ins;
		( p_ins.push_back( &ins ), ... ); // Fold expression
		return combine_channels( p_ins );
		}



	//============================================================================================================================================================
	// Information
	//============================================================================================================================================================

	/** Find the total energy in each channel of the input. Energy is the sum of the square of each sample. 
	 */
	std::vector<float> get_total_energy(
		) const;

	/** Find the difference in energy between two Audios. This is useful for unit testing algorithms.
	 */
	std::vector<float> get_energy_difference( 
		const Audio & other 
		) const;

	/** Try to find the wavelength of the input over time. Selects to minimize differences per repitition.
	 *	\param channel The channel to process.
	 *	\param start Frame to analyze.
	 *	\param absolute_cutoff Minimum the local wave difference minimum must achieve to be considered a potential wavelength.
	 *	\param window_size Input window size.
	 */
	float get_local_wavelength( 
		Channel channel,
		Frame start, 
		Frame window_size = 2048,
		float absolute_cutoff = 0.2f, 
		Frame minimum_wavelength = 10
		) const;

	/** Try to find the wavelength of the input over time. This is estimated using parabolic interpolation, so it may not be an integer.
	 *	\param channel The channel to process.
	 *	\param start Start frame.
	 *	\param end End frame.
	 *	\param absolute_cutoff Minimum the local wave difference minimum must achieve to be considered a potential wavelength.
	 *	\param window_size The amount of input data used in each analysis.
	 *	\param hop The hop size in frames per anlysis.
	 */
	std::vector<float> get_local_wavelengths( 
		Channel channel, 
		Frame start = 0, 
		Frame end = -1, 
		Frame window_size = 2048, 
		Frame hop = 128, 
		float absolute_cutoff = 0.2f,
		Frame minimum_wavelength = 10,
		flan_CANCEL_ARG 
		) const;

	/** Find the average of local wavelengths over time.
	 *	\param local_wavelengths Wavelengths to average. A -1 wavelength indicates no detectable wavelength. 
	 *	\param min_active_ratio This proportion of the inputs must have a detectable wavelength to continue computation. Otherwise, -1 is returned.
	 *	\param max_length_sigma If the standard deviation of the wavelengths is above this threshold, -1 is returned. Inputting -1 will disable sigma thresholding.
	 */
	float get_average_wavelength( 
		const std::vector<float> & local_wavelengths, 
		float min_active_ratio = 0, 
		float max_length_sigma = -1 
		) const; // Same as get_average_wavelengths

	/** Find the average of local wavelengths over time. This is a utility for calling get_local_wavelength and passing it to the other get_average_wavelength.
	 *	\param channel The channel to process.
	 *	\param min_active_ratio This proportion of the inputs must have a detectable wavelength to continue computation. Otherwise, -1 is returned.
	 *	\param max_length_sigma If the standard deviation of the wavelengths is above this threshold, -1 is returned. Inputting -1 will disable sigma thresholding.
	 *	\param start Start frame.
	 *	\param end End frame.
	 *	\param window_size The amount of input data used in each analysis.
	 *	\param hop The hop size in frames per anlysis.
	 */
	float get_average_wavelength( 
		Channel, 
		float min_active_ratio = 0, 
		float max_length_sigma = -1, 
		Frame start = 0, 
		Frame end = -1, 
		Frame window_size = 2048, 
		Frame hop = 128, 
		flan_CANCEL_ARG 
		) const; // Averages get_local_wavelengths
	
	/** Try to find the wavelength of the input over time. Selects highest to minimize octave error.
	 *	\param channel The channel to process.
	 *	\param start Frame to analyze.
	 *	\param window_size Input window size.
	 */
	Frequency get_local_frequency( 
		Channel channel, 
		Frame start = 0, 
		Frame window_size = 2048, 
		flan_CANCEL_ARG 
		) const;

	/** Utility for calling get_local_frequency at a set frame hop interval from start to end.
	 */
	std::vector<Frequency> get_local_frequencies( 
		Channel channel, 
		Frame start = 0, 
		Frame end = -1, 
		Frame window_size = 2048, 
		Frame hop = 128, 
		flan_CANCEL_ARG 
		) const;

	Function<Second, Amplitude> get_amplitude_envelope(
		Second window_width = 0.1f
		) const;

	Function<Second, Amplitude> get_frequency_envelope( 
		) const;



	//============================================================================================================================================================
	// Temporal
	//============================================================================================================================================================

	Audio modify_boundaries( 
		Second start, 
		Second end 
		) const;

	Audio remove_edge_silence( 
		Amplitude non_silent_level,
		Second fade_time = 0
		) const;

	Audio remove_silence(
		Amplitude non_silent_level,
		Second minimum_gap,
		Second fade_in_time = 0
		) const;

	/** This reverses the Audio.
	 */
	Audio reverse(
		) const;

	/** This returns a peice of the original Audio that lied between start_time and end_time.
	 *	\param start Start time for cut.
	 *	\param end End time for cut.
	 *	\param start_fade The cut audio start is faded for this amount of time.
	 *	\param end_fade The cut audio end is faded for this amount of time.
	 */
	Audio cut( 
		Second start, 
		Second end, 
		Second start_fade = 0, 
		Second end_fade = 0 
		) const;

	/** This returns a peice of the original Audio that lied between start and end.
	 *	\param start Start frame for cut.
	 *	\param end End frame for cut.
	 *	\param start_fade The cut audio start is faded for this amount of framea.
	 *	\param end_fade The cut audio end is faded for this amount of frames.
	 */
	Audio cut_frames( 
		Frame start, 
		Frame end, 
		Frame start_fade = 0, 
		Frame end_fade = 0 
		) const;

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
	 *	\param qual
	 */
	Audio repitch( 
		const Function<Second, float> & factor, 
		Second granularity = .001f, 
		WDLResampleType quality = WDLResampleType::Sinc 
		) const;

	/** This repeats the input n times.
	 *	\param num_iterations The number of desired iterations
	 *	\param mod This allows user-defined processing of each iteration.
	 *	\param feedback When enabled, the previous iteration is sent to mod, rather than the original Audio.
	 */
	Audio iterate( 
		int32_t num_iterations, 
		Second crossfade_time = 0,
		const AudioMod & mod = AudioMod(), 
		bool feedback = false 
		) const;

	/** Generates a volume-decaying iteration of the input. If a mod is supplied, it will 
	 *	be fed the previous iteration's ouput, rather than the original Audio.
	 *	\param length The time at which delay events stop generating.
	 *	\param delay_time The amount of time between delays.
	 *	\param decay Each delay event will have its gain scaled by the decay at the event time.
	 *	\param mod This allows user-defined processing of each delay. Delay line feedback will have this repeatedly applied.
	 */
	Audio delay( 
		Second length, 
		const Function<Second, Second> & delay_time, 
		const Function<Second, float> & decay = .5, 
		const AudioMod & mod = AudioMod() 
		) const;

	// Audio stereo_delay(
	// 	Second length,
	// 	const Function<Second, Second> & l_time,
	// 	const Function<Second, Second> & r_time,
	// 	const Function<Second, Amplitude> & decay
	// 	) const;

	std::vector<Audio> split_at_times(
		std::vector<Second> split_times, // Pass by value is intentional
		Second fade = 0
		) const; 

	std::vector<Audio> split_with_lengths(
		std::vector<Second> split_lengths, // Pass by value is intentional
		Second fade = 0
		) const;

	/** This slices the input into a vector of chunks, each slice_length seconds long.
	 *	The final chunk will be zero padded to meet the slice_length.
	 *  \param slice_length The length of each output.
	 *  \param fade The start and end fade time of each slice.
	 */
	std::vector<Audio> split_with_equal_lengths( 
		Second slice_length, 
		Second fade = 0
		) const;

	/** This performs split_with_equal_lengths, randomizes the order of the chopped pieces, and joins, crossfading.
	 *  \param slice_length The length of each output.
	 *  \param fade The start and end fade time of each slice.
	 */
	// Audio rearrange( 
	// 	Second slice_length, 
	// 	Second fade = 0
	// 	) const;

	Audio random_chunks(
		Second length,
		const Function<Second, Second> & chunk_length,
		const Function<Second, Second> & fade = 0,
		const AudioMod & mod = AudioMod()
		) const;

	//============================================================================================================================================================
	// Volume
	//============================================================================================================================================================

	/** Volume scaling, output( t ) = input( t ) * gain( t ).
	 *	\param gain Describes how volume should be scaled as a function of time.
	 */
	Audio modify_volume( 
		const Function<Second, float> & gain 
		) const;

	Audio ring_modulate( 
		const Audio & other 
		) const;

	Audio& modify_volume_in_place( 
		const Function<Second, float> & gain
		);		

	/** Volume setting. This normalizes (divides every sample by the largest sample value), and then scales by level.
	 *	\param level Describes how volume should be scaled as a function of time.
	 */
	Audio set_volume( 
		const Function<Second, Amplitude> & level 
		) const;

	Audio& set_volume_in_place(
		const Function<Second, Amplitude> & level
		);

	/** This adds a fade to the ends of the input Audio.
	 *	\param start Length of the start fade.
	 *	\param end Length of the end fade.
	 *	\param interp The fading curve. Square root should be used for constant-power crossfading.
	 */
	Audio fade( 
		Second start = 16.0f/48000.0f, 
		Second end = 16.0f/48000.0f, 
		Interpolator interp = Interpolators::sqrt 
		) const;

	Audio& fade_in_place(
	 	Second start = 16.0f/48000.0f, 
		Second end = 16.0f/48000.0f, 
		Interpolator interp = Interpolators::sqrt 
		);

	/** This adds a fade to the ends of the input Audio.
	 *	\param start Length of the start fade.
	 *	\param end Length of the end fade.
	 *	\param interp The fading curve. Square root should be used for constant-power crossfading.
	 */
	Audio fade_frames( 
		Frame start = 16, 
		Frame end = 16, 
		Interpolator interp = Interpolators::sqrt 
		) const;

	Audio& fade_frames_in_place( 
		Frame start = 16, 
		Frame end = 16, 
		Interpolator interp = Interpolators::sqrt 
		);

	/** This phase inverts the input. Every output sample is assigned the negative of the input sample at that time.
	 */
	Audio invert_phase(
		) const;

	/** This applies the shaper as a function to each sample in the input.
	 *	\param shaper Each sample in the input is passed through this. 
	 *		Samples will be values on [-1,1] under normal circumstances.
	 *		For example, the function y = x would have no effect as a shaper.
	 */
	Audio waveshape( 
		const Function< std::pair<Second, Sample>, Sample > & shaper,
		uint16_t oversample_factor = 4
		) const;

	/** This is meant to be a black box process for adding a moisture effect to bass signals.
	 *  \param amount How much of the effect to add.
	 *  \param frequency Base effect frequency when skew is 1.
	 *  \param skew Magic ingredient.
	 */
	Audio add_moisture(
		const Function<Second, Amplitude> & amount = .5f,
		const Function<Second, Frequency> & frequency = 96,
		const Function<Second, float> & skew = 4,
		const Function<Second, Amplitude> & waveform = waveforms::sine
		) const;

	/** This is a dynamic range compressor.
	 */
	Audio compress( 
		const Function<Second, Decibel> & threshold, 
		const Function<Second, float> & compression_ratio = 3.0f, 
		const Function<Second, Second> & attack = 5.0f / 1000.0f, 
		const Function<Second, Second> & release = 100.0f / 1000.0f, 
		const Function<Second, Decibel> & knee_width = Decibel( 0 ), 
		const Audio * sidechain_source = nullptr 
		) const;

	Audio apply_adsr_envelope(
		Second attack_time,
		Second decay_time,
		Second sustain_time,
		Second release_time,
		Amplitude sustain_level,
		float attack_exponent = 1,
		float decay_exponent = 1,
		float release_exponent = 1
		) const;


	// Helper for common adsr usage
	Audio apply_ar_envelope(
		Second attack_time,
		Second release_time,
		float attack_exponent = 1,
		float release_exponent = 1
		) const;

	//============================================================================================================================================================
	// Spatial
	//============================================================================================================================================================

	/// This can't be implemented until I have a setup with more than two speakers to test it on
	// /** This spatially repositions the input. This is a general purpose panning algorithm for handling any number of speakers. 
	//  *	The speakers are assumed to be in a plane. Units are given in meters for haas effect purposes.
	//  *	\param pan_position The plane position of the input over time. 
	//  *  \param listenerPosition The plane position of the listener over time.
	//  *	\param speakerPositions The plane positions of each speaker.
	//  *  \param hassEffect Decides if the haas effect should be applied. 
	//  */
	// Audio pan( Function<Second,vec2> pan_position, Function<Second, vec2> listenerPosition, std::vector<vec2> speakerPositions, bool haasEffect = false ) const;

	/// @brief This uses a bunch of psychoacoustic techniques to give a 3d spatial postion to a mono input.
	/// @param position Where the sound lives over time.
	/// @return 
	Audio stereo_spatialize( 
		const Function<Second, vec2> & position,
		Meter head_width = 0.18f
		) const;

	/** This spatially repositions the input for mono and stereo inputs. Stereo speakers are assumed form an equilateral trangle with side
	 *	length 1 meter with the listener. The input position is assumed to sit halfway between the speakers.
	 *  If the number of channels is neither 1 nor 2 a null Audio is returned.
	 *	\param pan_position The input spatial position from -1 (left) to 1 (right) over time.
	 */
	Audio pan( 
		const Function<Second, float> & pan_position 
		) const;

	Audio& pan_in_place( 
		const Function<Second, float> & pan_position 
		);

	/** This redistributes energy between the mid and side signals
	 *	\param widen_amount This should return a value on [-1,1], representing movement from mid to side.
	 */
	Audio widen( 
		const Function<Second, float> & widen_amount 
		) const;



	//============================================================================================================================================================
	// Filters
	//============================================================================================================================================================

	/*
	There are an essentially unlimited number of filter types and use cases, so I have made no attempt to cover all of them.
	If you are reading this and your experience with filters is purely musical, you will find 2-pole filters to be
	the familiar resonant filter type found in most DAW programs. The "bell" shaped filter found in most EQs is here called
	a "band shelving filter".

	Everything to do with filters here is based on the book "VA Filter Design" (second edition).
	https://ia601900.us.archive.org/5/items/the-art-of-va-filter-design-rev.-2.1.2/VAFilterDesign_2.1.2.pdf#chapter.10
	This textbook can be very rough around the edges, but most questions on why I am building filters a certain way can
	be answered by searching it.
	*/

	/* 1 Pole ============================ */

	Audio filter_1pole_lowpass(
		const Function<Second, Frequency> & cutoff,
		uint16_t order = 1
		) const;

	Audio filter_1pole_highpass(
		const Function<Second, Frequency> & cutoff,
		uint16_t order = 1
		) const;

	std::vector<Audio> filter_1pole_split(
		const Function<Second, Frequency> & cutoff,
		uint16_t order = 1
		) const;

	Audio filter_1pole_lowshelf(
		const Function<Second, Frequency> & cutoff,
		const Function<Second, Decibel> & gain,
		uint16_t order = 1
		) const;

	Audio filter_1pole_highshelf(
		const Function<Second, Frequency> & cutoff,
		const Function<Second, Decibel> & gain,
		uint16_t order = 1
		) const;

	/** This applies the same 1-pole low pass filter to the input n times. 
		It is mainly a tool used for modeling atmospheric scattering in spacialization methods.
	*/
	Audio filter_1pole_repeat_low(
		const Function<Second, Frequency> & cutoff,
		const uint16_t repeats
		) const;

	Audio filter_1pole_repeat_high(
		const Function<Second, Frequency> & cutoff,
		const uint16_t repeats
		) const;

	/* 2 Pole ============================ */

	Audio filter_2pole_lowpass(
		const Function<Second, Frequency> & cutoff,
		const Function<Second, float> & damping,
		uint16_t order = 1
		) const;

	Audio filter_2pole_bandpass(
		const Function<Second, Frequency> & cutoff,
		const Function<Second, float> & damping,
		uint16_t order = 1
		) const;

	Audio filter_2pole_highpass(
		const Function<Second, Frequency> & cutoff,
		const Function<Second, float> & damping,
		uint16_t order = 1
		) const;

	std::vector<Audio> filter_2pole_split(
		const Function<Second, Frequency> & cutoff,
		const Function<Second, float> & damping,
		uint16_t order = 1
		) const;

	Audio filter_2pole_notch(
		const Function<Second, Frequency> & cutoff,
		const Function<Second, float> & damping,
		uint16_t order = 1
		) const;

	Audio filter_2pole_lowshelf(
		const Function<Second, Frequency> & cutoff,
		const Function<Second, float> & damping,
		const Function<Second, Decibel> & gain,
		uint16_t order = 1
		) const;

	Audio filter_2pole_bandshelf(
		const Function<Second, Frequency> & cutoff,
		const Function<Second, float> & damping,
		const Function<Second, Decibel> & gain,
		uint16_t order = 1
		) const;

	Audio filter_2pole_highshelf(
		const Function<Second, Frequency> & cutoff,
		const Function<Second, float> & damping,
		const Function<Second, Decibel> & gain,
		uint16_t order = 1
		) const;

	/* Other filters ==================== */

	Audio filter_1pole_multinotch(
		const Function<Second, Frequency> & cutoff,
		const Function<Second, float> & wet_dry = .5,
		uint16_t order = 1,
		bool invert = false
		) const;

	Audio filter_2pole_multinotch(
		const Function<Second, Frequency> & cutoff,
		const Function<Second, float> & damping,
		const Function<Second, float> & wet_dry = .5,
		uint16_t order = 1,
		bool invert = false
		) const;

	Audio filter_comb(
		const Function<Second, Frequency> & cutoff,
		const Function<Second, float> & feedback = 0,
		const Function<Second, float> & wet_dry = .5,
		bool invert = false
		) const;



	//============================================================================================================================================================
	// Combination
	//============================================================================================================================================================

	/** The main Audio mixing method. The ouput format will be that of the first input.
	 *	\param ins The Audio to mix.
	 *	\param start_times This gives the start time for each input Audio. The default will start all inputs at time 0.
	 *	\param balances The volume scaling to apply to each input Audio. The default will apply no scaling.
	 *		These are relative to global time 0 rather than the start time of each input.
	 */
	static Audio mix( 
		std::vector<const Audio *> ins, 
		std::vector<Second> start_times = {}, 
		std::vector<Amplitude> gains = {} 
		);

	static Audio mix( 
		const std::vector<Audio> & ins, 
		std::vector<Second> start_times = {}, 
		std::vector<Amplitude> gains = {} 
		);

	template<typename ... Ts, typename = AllAudios<Ts...> >
	static Audio mix( 
		const Ts & ... ins 
		)
		{
		std::vector<const Audio *> p_ins;
		( p_ins.push_back( &ins ), ... ); // Fold expression
		return mix( p_ins );
		}

	static Audio mix_variable_gain( 
		std::vector<const Audio *> ins, 
		std::vector<Second> start_times, 
		const std::vector<const Function<Second, Amplitude> *> & gains 
		);

	static Audio mix_variable_gain( 
		const std::vector<Audio> & ins, 
		const std::vector<Second> & start_times, 
		const std::vector<Function<Second, Amplitude>> & gains 
		);

	Audio& mix_in_place( 
		const Audio & other, 
		Second other_start_time = 0, 
		const Function<Second, Amplitude> & other_amplitude = 1.0f 
		);

	static Audio join( 
		const std::vector<const Audio *> & ins, 
		const std::vector<Second> & offsets
		);

	/** This joins all input Audio tip to tail in the order they were passed.
	 *	\param ins The Audio to join.
	 *	\param offset 
	 */
	static Audio join( 
		const std::vector<const Audio *> & ins, 
		Second offset = 0 
		);

	static Audio join( 
		const std::vector<Audio> & ins, 
		Second offset = 0 
		);

	template<typename ... Ts, typename = AllAudios<Ts...> >
	static Audio join( 
		const Ts & ... ins 
		)
		{
		std::vector<const Audio *> p_ins;
		( p_ins.push_back( &ins ), ... ); // Fold expression
		return Audio::join( p_ins );
		}

	/** At each point in time, selection decides which of the input Audio streams is playing.
	 *	Non-integer selections will mix appropriately scaled copies of the surrounding integer inputs.
	 *	\param ins The Audio to mix.
	 *	\param selection Which audio to play.
	 *	\param start_times Second offsets for each input.
	 */
	static Audio select( 
		const std::vector<const Audio *> & ins, 
		const Function<Second, float> & selection, 
		const std::vector<Second> & start_times = {} 
		);

	static Audio select( 
		const std::vector<Audio> & ins, 
		const Function<Second, float> & selection, 
		std::vector<Second> start_times = {} 
		);

	template<typename ... Ts, typename = AllAudios<Ts...> >
	static Audio select( 
		const Function<Second, float> & selection, 
		const Ts & ... ins 
		)
		{
		std::vector<const Audio *> p_ins;
		( p_ins.push_back( &ins ), ... ); // Fold expression
		return select( p_ins, selection );
		}

	/** This is normal convolution between Audio. If the number of channels are mismatched, the ir channels will be used cyclically.
	 *	\param ir An impulse response.
	 */
	Audio convolve( 
		const Audio & ir,
		bool normalize = true
		) const;



	//============================================================================================================================================================
	// Synthesis
	//============================================================================================================================================================

	/** Generate an Audio from a waveform function. This can be accomplished in a more general setting using the function ctor of flan::Wavetable.
	*  \param waveform This is evaluated from 0 to 1. This portion is repeated to create the Audio.
	*  \param length The length of the output.
	*  \param freq The frequency of the output.
	*  \param sample_rate The sample rate of the output.
	*  \param oversample Synthesized audio must be generated at a sample rate higher than the desired sample rate to avoid aliasing.
	*      This describes how many samples should be used in the synthesis per sample in the output.
	*/
	static Audio synthesize_waveform(  
		const Function<Second, Amplitude> & waveform, 
		Second length, 
		const Function<Second, Frequency> & freq, 
		FrameRate sample_rate = 48000, 
		size_t oversample = 16 
		);

	static Audio synthesize_white_noise(
		Second length,
		FrameRate sample_rate = 48000, 
		size_t oversample = 16 
		);

	static Audio synthesize_pink_noise(
		Second length,
		FrameRate sample_rate = 48000,
		int num_rows = 128
		);

	/** This process generates an Audio given information about the spectral dispersion and amplitude of each harmonic.
	*  \param length The length of the output.
	*  \param freq The frequency of the output.
	*  \param spread The standard deviation of each harmonic peak if using a normal distribution. Scales as a standard deviation would for exotic distributions.
	*  \param harmonic_scale A scaling factor on each harmonic.
	*  \param peak_distribution The spectral energy distribution within each harmonic.
	*  \param fundamental_power Raising two to this power gives the spectral space between harmonics.
	*  \param spectrum_size_power Raising two to this power gives the total spectrum size. Maximum of 25.
	*  \param granularity The frequency at which input functions are sampled.
	*/
	static Audio synthesize_spectrum(
		Second length, 
        const Function<Second, Frequency> & freq,
		const Function<Harmonic, Frequency> & spread = []( Harmonic h ){ return h; },
		const Function<Harmonic, Amplitude> & harmonic_scale = []( Harmonic h ){ return 1.0f/std::sqrt(h); },
		const Function<Frequency, Amplitude> & peak_distribution = []( Frequency x )
			{ 
			return 1.0f / std::sqrt(pi2) * std::exp( -0.5f * std::pow( x, 2.0f ) );
			},
		int fundamental_power = 8,
		int spectrum_size_power = 20,
		Channel num_channels = 2,
        Second granularity = 0.001f
		);

	static Audio synthesize_impulse(
		Frequency base_freq,
		Harmonic num_harmonics = std::numeric_limits<Harmonic>::max(),
		float chroma = 1.0f,
		FrameRate sample_rate = 48000.0f
		);

	// Grain Controllers ===================================================================================================================

	static Audio synthesize_grains( 
		Second length,
		const Function<Second, float> & grains_per_second,
		const Function<Second, float> & time_scatter,
		const Function<Second, Audio> & grain_source,
		FrameRate sample_rate = 48000
		);

	/* Granular synthesis using a single input Audio as the grain source, optionally passed through a mod function.
	 */
	Audio texture(
		Second length, 
		const Function<Second, float> & grains_per_second, 
		const Function<Second, float> & time_scatter, 
		const AudioMod & mod = AudioMod(),
		bool mod_feedback = false
		) const;

	// Grain Compositions ===================================================================================================================

	/* Trainlet synthesis as seen in the book "Microsound".
	 * 	\param grains_per_second See Audio::synthesize_grains.
	 * 	\param time_scatter See Audio::synthesize_grains.
	 * 	\param position Spatial plane position to spacialize trainlet at.
	 * 	\param trainlet_gain_envelope The amplitude over time of individual trainlets 
	 * 	\param impulse_freq The rate at which impulses occur within an individual trainlet. 
	 * 	\param trainlet_length 
	 * 	\param num_harmonics How many harmonics should be used to synthesize an impulse.
	 * 	\param chroma A geometric scaling factor on the harmonics within an impulse.
	 * 	\param impulse_harmonic_frequency The frequency of the harmonics used in impulse generation.
	 * 	\param sample_rate
	 */
	static Audio synthesize_trainlets( 
		Second length,
		const Function<Second, float> & grains_per_second,
		const Function<Second, float> & time_scatter,
		const Function<Second, vec2> & position,
		const Function<Second, Amplitude> & trainlet_gain_envelope,
		const Function<Second, Frequency> & impulse_freq,
		const Function<Second, Second> & trainlet_length,
		const Function<Second, Harmonic> & num_harmonics = std::numeric_limits<Harmonic>::max(),
		const Function<Second, float> & chroma = 1,
		const Function<Second, Frequency> & impulse_harmonic_frequency = 32,
		FrameRate sample_rate = 48000
		);

	/** Granular synthesis using sections of a single input as the grain source.
	 *  \param length Grains will generate until this time.
	 *  \param grains_per_second The mean number of grains per second.
	 *  \param scatter The standard deviation in grains. Due to being in grains, a higher 
	 *		grains_per_second will cause scatter to have less of an effect.
	 *  \param selection The input time from which grains will be read.
	 *  \param grain_length The length of grains.
	 *  \param fade_time The start and end fade time for each grain. Fades use a sqrt curve.
	 *  \param mod Each grain is fed into this along with the event time.
	 */
	Audio granulate( 
		Second length, 
		const Function<Second, float> & grains_per_second, 
		const Function<Second, float> & time_scatter, 
		const Function<Second, Second> & time_selection, 
		const Function<Second, Second> & grain_length,
		const Function<Second, Second> & fade_time = 0,
		const AudioMod & mod = AudioMod()
		) const;

	Audio psola(
		Second length, 
		const Function<Second, Second> & time_selection,
		const AudioMod & mod = AudioMod()
		) const;

	// static Audio synthesize_pulsars(
	// 	Second length,
	// 	const Function<Second, Frequency> & pulse_frequency, 
	// 	const Function<Second, Amplitude> & waveform, 
	// 	const Function<Second, Frequency> & waveform_frequency,
	// 	const Function<float, Amplitude> & pulsaret_envelope,
	// 	);

private:
	static std::vector<Audio> match_sample_rates_or_return_null( const std::vector<const Audio *> & ins );

	};

} // End namespace flan
