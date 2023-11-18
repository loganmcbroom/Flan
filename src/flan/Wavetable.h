#pragma once

#include "flan/Audio/Audio.h"

namespace flan {

/** -1 will mean "Apply to all", otherwise should be a value on [0,1] */
using Waveform = float;

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

    /** Main Audio wavetable constructor.
     * \param source The Audio source.
     * \param snap_mode This decides how waveforms starts and ends are nudged to avoid edge discontinuities. See Waveform::SnapMode.
     * \param pitch_mode This decides how waveforms lengths are estimated. See Waveform::PitchMode.
     * \param snap_ratio The proportion of the expected frame jump to search for frame snapping.
     * \param fixed_frame_size When pitch_mode is None, this is used in place of pitch deduction.
     */
    Wavetable( 
        const Audio & source, 
        SnapMode snap_mode = SnapMode::Zero, 
        PitchMode pitch_mode = PitchMode::Local, 
        float snap_ratio = .3, 
        Frame fixed_frame_size = 256, 
        flan_CANCEL_ARG );

    /** Functional wavetable constructor.
     * \param f The function to sample. Each wave will be sampled on the interval [k,k+1), starting at 0.
     * \param num_waves The number of waveforms to sample.
     */
    Wavetable( 
        const Function<Second, Amplitude> & f, 
        int num_waves, 
        flan_CANCEL_ARG );

    // /** This constructs a Wavetable using specific waveform start points. This is normally used only by the main constructor which calls get_waveform_starts. 
    //  *  If get_waveform_starts is having a hard time, you can manually set the waveform start frames and use this.
    //  * \param source The audio source.
    //  * \param waveform_starts The start frames of each waveform. The waveforms are assumed to be contiguous.
    //  */
    // Wavetable( 
    //     const Audio & source, 
    //     const std::vector<std::vector<Frame>> & waveform_starts, 
    //     flan_CANCEL_ARG );

    /** Generate Audio.
     * \param length The output length.
     * \param freq The output frequency over time.
     * \param index The waveform to read from. Values from zero to one will read from the start to the end of the table.
     * \param granularity The time between input function evaluations.
     */
    Audio synthesize( 
        Second length, 
        const Function<Second, Frequency> & freq, 
        const Function<Second, Waveform> & index,
        bool smooth = true, 
        Second granularity = 0.001f, 
        flan_CANCEL_ARG 
        ) const;

    Graph graph_waveform( int waveform_index ) const;

    Graph graph_waveform_range( Channel channel, int start, int end ) const;

    int get_num_waveforms( Channel ) const;

    /** Fade the edges of target waveforms towards zero.
     * \param fade_frames The fade length.
     */
    void add_fades_in_place( Frame fade_frames = 32 );

    /** Fade the edges of target waveforms towards the halfway point between the edge heights.
     * \param fade_frames The fade length.
     */
    void remove_jumps_in_place( Frame fade_frames = 32 );

    /** Remove the DC offset from target waveforms.
     */
    void remove_dc_in_place();

    void normalize_in_place();

private:
    Wavetable();

    Sample * get_waveform( int waveform_index, Channel channel );
    const Sample * get_waveform( int waveform_index, Channel channel ) const;

    float ratio_to_table_index( float r, Channel ) const;

    bool is_null() const;

    const Frame wavelength;
    const fFrame num_source_frames;
    // The outer vector is for the channel you want to access.
    // The inner vector keeps what source frame each waveform was drawn from
    std::vector<std::vector<Frame>> waveform_starts;

    Audio table;
};

};
