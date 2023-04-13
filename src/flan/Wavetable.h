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

    static std::vector<Frame> getWaveformStarts( const Audio & source, Wavetable::SnapMode snapMode, Wavetable::PitchMode pitchMode, 
        float snapRatio, Frame fixedFrame, flan_CANCEL_ARG );

    /** Main Audio wavetable constructor.
     * \param source The Audio source.
     * \param snapMode This decides how waveforms starts and ends are nudged to avoid edge discontinuities. See Waveform::SnapMode.
     * \param pitchMode This decides how waveforms lengths are estimated. See Waveform::PitchMode.
     * \param snapRatio The proportion of the expected frame jump to search for frame snapping.
     * \param fixedFrameSize When pitchMode is None, this is used in place of pitch deduction.
     */
    Wavetable( const Audio & source, SnapMode snapMode = SnapMode::Zero, PitchMode pitchMode = PitchMode::Local, float snapRatio = .3, 
        Frame fixedFrameSize = 256, flan_CANCEL_ARG );

    /** Functional wavetable constructor.
     * \param f The function to sample. Each wave will be sampled on the interval [k,k+1), starting at 0.
     * \param numWaves The number of waveforms to sample.
     */
    Wavetable( Func1x1 f, int numWaves, flan_CANCEL_ARG );

    /** This constructs a Wavetable using specific waveform start points. This is normally used only by the main constructor which calls getWaveformStarts. 
     *  If getWaveformStarts is having a hard time, you can manually set the waveform start frames and use this.
     * \param source The audio source.
     * \param waveformStarts The start frames of each waveform. The waveforms are assumed to be contiguous.
     */
    Wavetable( const Audio & source, const std::vector<Frame> & waveformStarts, flan_CANCEL_ARG );

    /** Generate Audio.
     * \param length The output length.
     * \param freq The output frequency over time.
     * \param index The waveform to read from. Choosing not-integer indices will read a linear interpolation of the surrounding waveforms.
     * \param granularity The time between input function evaluations.
     */
    Audio synthesize( Time length, Func1x1 freq, Func1x1 index, Time granularity = 0.001f, flan_CANCEL_ARG ) const;

    /** Fade the edges of target waveforms towards zero.
     * \param fadeFrames The fade length.
     * \param target The waveform to modify. Use -1 to target all.
     */
    void addFades( Frame fadeFrames = 32, Waveform target = -1 );

    /** Fade the edges of target waveforms towards the halfway point between the edge heights.
     * \param fadeFrames The fade length.
     * \param target The waveform to modify. Use -1 to target all.
     */
    void removeJumps( Frame fadeFrames = 32, Waveform target = -1 );

    /** remove the DC offset from target waveforms.
     * \param target The waveform to modify. Use -1 to target all.
     */
    void removeDC( Waveform target = -1 );

private:
    Wavetable();

    Sample * getWaveform( Waveform waveform );

    const Frame wavelength;
    Audio table;
    const int numWaveforms;
};

};
