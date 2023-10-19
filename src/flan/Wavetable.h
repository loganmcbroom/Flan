#pragma once

#include "flan/Audio/Audio.h"

namespace flan {

/** -1 will imply "Apply to all" */
using Waveform = int;

/* Functions for working with wavetable synthesis
 */
class Wavetable
{
public:

    enum class SnapMode
        {
        None, /** Make no attempt to snap frames. */
        Zero, /** Waveforms will try to begin and end at zero-crossings. */
        Level, /** Waveforms will try to begin and end at the same level. */
        };

    enum class PitchMode
        {
        None, /** The input pitch isn't estimated. The supplied frame size is always used. */
        Local, /** The input pitch is estimated as a function of time. */
        Global, /** The average pitch of the input is estimated. */
        };

    static std::vector<Frame> get_waveform_starts( const Audio & source, Wavetable::SnapMode snap_mode, Wavetable::PitchMode pitch_mode, 
        float snap_ratio, Frame fixed_frame, flan_CANCEL_ARG );

    /** Main Audio wavetable constructor.
     * \param source The Audio source.
     * \param snap_mode This decides how waveforms starts and ends are nudged to avoid edge discontinuities. See Waveform::SnapMode.
     * \param pitch_mode This decides how waveforms lengths are estimated. See Waveform::PitchMode.
     * \param snap_ratio The proportion of the expected frame jump to search for frame snapping.
     * \param fixed_frameSize When pitch_mode is None, this is used in place of pitch deduction.
     */
    Wavetable( const Audio & source, SnapMode snap_mode = SnapMode::Zero, PitchMode pitch_mode = PitchMode::Local, float snap_ratio = .3, 
        Frame fixed_frameSize = 256, flan_CANCEL_ARG );

    /** Functional wavetable constructor.
     * \param f The function to sample. Each wave will be sampled on the interval [k,k+1), starting at 0.
     * \param numWaves The number of waveforms to sample.
     */
    Wavetable( Function<Second, Amplitude> f, int numWaves, flan_CANCEL_ARG );

    /** This constructs a Wavetable using specific waveform start points. This is normally used only by the main constructor which calls get_waveform_starts. 
     *  If get_waveform_starts is having a hard time, you can manually set the waveform start frames and use this.
     * \param source The audio source.
     * \param waveform_starts The start frames of each waveform. The waveforms are assumed to be contiguous.
     */
    Wavetable( const Audio & source, const std::vector<Frame> & waveform_starts, flan_CANCEL_ARG );

    /** Generate Audio.
     * \param length The output length.
     * \param freq The output frequency over time.
     * \param index The waveform to read from. Choosing not-integer indices will read a linear interpolation of the surrounding waveforms.
     * \param granularity The time between input function evaluations.
     */
    Audio synthesize( Second length, Function<Second, Frequency> freq, Function<Second, float> index, Second granularity = 0.001f, flan_CANCEL_ARG ) const;

    /** Fade the edges of target waveforms towards zero.
     * \param fade_frames The fade length.
     * \param target The waveform to modify. Use -1 to target all.
     */
    void add_fades_in_place( Frame fade_frames = 32, Waveform target = -1 );

    /** Fade the edges of target waveforms towards the halfway point between the edge heights.
     * \param fade_frames The fade length.
     * \param target The waveform to modify. Use -1 to target all.
     */
    void remove_jumps_in_place( Frame fade_frames = 32, Waveform target = -1 );

    /** remove the DC offset from target waveforms.
     * \param target The waveform to modify. Use -1 to target all.
     */
    void remove_dc_in_place( Waveform target = -1 );

private:
    Wavetable();

    Sample * getWaveform( Waveform waveform );

    const Frame wavelength;
    Audio table;
    const int numWaveforms;
};

};
