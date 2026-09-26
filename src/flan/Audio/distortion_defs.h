#pragma once

#include <array>
#include <optional>

#include "filter_defs.h"

namespace flan {

struct Audio;

struct Shaper
    {
    virtual ~Shaper() = default;
    virtual void init( const Audio& ) {}
    virtual void reset() = 0; 
    virtual Sample process(Second t, Sample s) = 0; 
    };

struct FunctionShaper : public Shaper
    {
    Function<std::pair<Second, Sample>, Sample> function;

    FunctionShaper( const Function<std::pair<Second, Sample>, Sample>& function_ )
        : function( function_.copy() )
        {}
    void reset() override {}
    Sample process(Second t, Sample s) override { return function( std::pair( t, s ) ); }
    };

struct GainShaper : public Shaper
    {
    Function<Second, float> gain;

    GainShaper( const Function<Second, float>& gain ) : gain( gain.copy() ) {}
    void reset() override {}
    Sample process(Second t, Sample s) override;
    };

struct SoftClipShaper : public Shaper
    {
    Function<Second, float> drive;

    SoftClipShaper( const Function<Second, float>& drive ) : drive( drive.copy() ) {}
    void reset() override {}
    Sample process(Second t, Sample s) override;
    };

struct RectifierShaper : public Shaper
    {
    bool full_wave;

    RectifierShaper( bool full_wave = true ) : full_wave( full_wave ) {}
    void reset() override {}
    Sample process(Second t, Sample s) override;
    };

struct QuantizeShaper : public Shaper
    {
    Function<Second, Sample> step;

    QuantizeShaper( const Function<Second, Sample>& step ) : step( step.copy() ) {}
    void reset() override {}
    Sample process(Second t, Sample s) override;
    };

struct HardClipShaper : public Shaper
    {
    Function<Second, float> gain;

    HardClipShaper( const Function<Second, float>& gain ) : gain( gain.copy() ) {}
    void reset() override {}
    Sample process(Second t, Sample s) override;
    };

struct DCOffsetShaper : public Shaper
    {
    Function<Second, Sample> offset;

    DCOffsetShaper( const Function<Second, Sample>& offset ) : offset( offset.copy() ) {}
    void reset() override {}
    Sample process(Second t, Sample s) override;
    };

struct SineFoldShaper : public Shaper
    {
    Function<Second, float> drive;

    SineFoldShaper( const Function<Second, float>& drive ) : drive( drive.copy() ) {}
    void reset() override {}
    Sample process(Second t, Sample s) override;
    };

struct CubicShaper : public Shaper
    {
    Function<Second, float> drive;

    CubicShaper( const Function<Second, float>& drive ) : drive( drive.copy() ) {}
    void reset() override {}
    Sample process(Second t, Sample s) override;
    };

struct ExponentialShaper : public Shaper
    {
    Function<Second, float> drive;

    ExponentialShaper( const Function<Second, float>& drive ) : drive( drive.copy() ) {}
    void reset() override {}
    Sample process(Second t, Sample s) override;
    };

struct WavefolderShaper : public Shaper
    {
    Function<Second, float> drive;

    WavefolderShaper( const Function<Second, float>& drive ) : drive( drive.copy() ) {}
    void reset() override {}
    Sample process(Second t, Sample s) override;
    };

struct TremoloShaper : public Shaper
    {
    Function<Second, Frequency> frequency;
    Function<Second, float> depth;
    FrameRate sr;
    float phase = 0;

    TremoloShaper(
        const Function<Second, Frequency>& frequency,
        const Function<Second, float>& depth = 1.0f
        )
        : frequency( frequency.copy() ), depth( depth.copy() ) {}
    void init( const Audio& me ) override;
    void reset() override;
    Sample process(Second t, Sample s) override;
    };

struct DownsampleShaper : public Shaper
    {
    Function<Second, FrameRate> hold_rate;
    FrameRate sr;
    float phase = 0;
    Sample held_sample = 0;
    bool initialized = false;

    DownsampleShaper( const Function<Second, FrameRate>& hold_rate )
        : hold_rate( hold_rate.copy() ) {}
    void init( const Audio& me ) override;
    void reset() override;
    Sample process(Second t, Sample s) override;
    };

struct RingModShaper : public Shaper
    {
    Function<Second, Frequency> frequency;
    Function<Second, float> amount;
    FrameRate sr;
    float phase = 0;

    RingModShaper(
        const Function<Second, Frequency>& frequency,
        const Function<Second, float>& amount = 1.0f
        )
        : frequency( frequency.copy() )
        , amount( amount.copy() ) 
        {}
    void init( const Audio& me ) override;
    void reset() override;
    Sample process(Second t, Sample s) override;
    };

struct SlewShaper : public Shaper
    {
    Function<Second, Sample> rise_rate;
    Function<Second, Sample> fall_rate;
    FrameRate sr;
    Sample output = 0;

    SlewShaper(
        const Function<Second, Sample>& rise_rate,
        const Function<Second, Sample>& fall_rate
        )
        : rise_rate( rise_rate.copy() ), fall_rate( fall_rate.copy() ) {}
    void init( const Audio& me ) override;
    void reset() override;
    Sample process(Second t, Sample s) override;
    };

struct AllpassShaper : public Shaper
    {
    Function<Second, Second> delay_length;
    Function<Second, float> feedback;
    Second max_delay_time;
    std::vector<Sample> input_buffer;
    std::vector<Sample> output_buffer;
    Frame write_pos = 0;
    FrameRate sr;

    AllpassShaper(
        const Function<Second, Second>& delay_length,
        const Function<Second, float>& feedback,
        Second max_delay_time = 1.0f
        )
        : delay_length( delay_length.copy() )
        , feedback( feedback.copy() )
        , max_delay_time( max_delay_time ) {}
    void init( const Audio& me ) override;
    void reset() override;
    Sample process(Second t, Sample s) override;
    };

struct CombShaper : public Shaper
    {
    Function<Second, Frequency> frequency;
    Function<Second, float> feedback;
    Frequency min_frequency;
    std::vector<Sample> buffer;
    Frame write_pos = 0;
    FrameRate sr;

    CombShaper(
        const Function<Second, Frequency>& frequency,
        const Function<Second, float>& feedback,
        Frequency min_frequency = 1.0f
        )
        : frequency( frequency.copy() )
        , feedback( feedback.copy() )
        , min_frequency( min_frequency ) {}
    void init( const Audio& me ) override;
    void reset() override;
    Sample process(Second t, Sample s) override;
    };

struct EnvelopeFollowerShaper : public Shaper
    {
    Function<Second, float> amount;
    Function<Second, Second> attack;
    Function<Second, Second> release;
    FrameRate sr;
    Sample envelope = 0;

    EnvelopeFollowerShaper(
        const Function<Second, float>& amount,
        const Function<Second, Second>& attack = 0.001f,
        const Function<Second, Second>& release = 0.1f
        )
        : amount( amount.copy() )
        , attack( attack.copy() )
        , release( release.copy() ) {}
    void init( const Audio& me ) override;
    void reset() override;
    Sample process(Second t, Sample s) override;
    };

struct CompressorShaper : public Shaper
    {
    Function<Second, float> ratio;
    Function<Second, float> threshold;
    Function<Second, Second> attack;
    Function<Second, Second> release;
    FrameRate sr;
    Sample envelope = 0;

    CompressorShaper(
        const Function<Second, float>& ratio,
        const Function<Second, float>& threshold = 0.5f,
        const Function<Second, Second>& attack = 0.001f,
        const Function<Second, Second>& release = 0.1f
        )
        : ratio( ratio.copy() )
        , threshold( threshold.copy() )
        , attack( attack.copy() )
        , release( release.copy() ) {}
    void init( const Audio& me ) override;
    void reset() override;
    Sample process(Second t, Sample s) override;
    };

struct DiodeShaper : public Shaper
    {
    Function<Second, float> drive;
    Function<Second, float> asymmetry;

    DiodeShaper(
        const Function<Second, float>& drive,
        const Function<Second, float>& asymmetry = 0.0f
        )
        : drive( drive.copy() ), asymmetry( asymmetry.copy() ) {}
    void reset() override {}
    Sample process(Second t, Sample s) override;
    };

struct NoiseGateShaper : public Shaper
    {
    Function<Second, float> threshold;
    Function<Second, float> range;
    Function<Second, Second> attack;
    Function<Second, Second> release;
    FrameRate sr;
    Sample envelope = 0;

    NoiseGateShaper(
        const Function<Second, float>& threshold,
        const Function<Second, float>& range = 0.0f,
        const Function<Second, Second>& attack = 0.001f,
        const Function<Second, Second>& release = 0.05f
        )
        : threshold( threshold.copy() )
        , range( range.copy() )
        , attack( attack.copy() )
        , release( release.copy() ) {}
    void init( const Audio& me ) override;
    void reset() override;
    Sample process(Second t, Sample s) override;
    };

// struct AsymmetricClipShaper : public Shaper
//     {
//     Function<Second, float> positive_limit;
//     Function<Second, float> negative_limit;
//     Function<Second, float> gain;

//     AsymmetricClipShaper(
//         const Function<Second, float>& positive_limit,
//         const Function<Second, float>& negative_limit,
//         const Function<Second, float>& gain = 1.0f
//         )
//         : positive_limit( positive_limit.copy() )
//         , negative_limit( negative_limit.copy() )
//         , gain( gain.copy() ) {}
//     void reset() override {}
//     Sample process(Second t, Sample s) override;
//     };

struct Shaper_1Pole 
    {
    int order;
    Function<Second, Frequency> cutoff;
    bool lowpass;
	std::vector<Pole> poles;

    FrameRate sr;
	std::optional<Filter_1Pole> filter_1pole;
	std::vector<Filter_2Pole> filter_2poles; 

    Shaper_1Pole( int order, const Function<Second, Frequency>& cutoff, bool lowpass );
    void init( const Audio& me );
    void reset();
    Sample process(Second frame, Sample sample);
    };

struct LowpassShaper : public Shaper
    {
    Shaper_1Pole shaper;
    LowpassShaper( int order, const Function<Second, Frequency>& cutoff ) : shaper( order, cutoff, true ) {}
    void init( const Audio& me ) override { shaper.init( me ); }
    void reset() override { shaper.reset(); }
    Sample process(Second t, Sample s) override { return shaper.process( t, s ); }
    };

struct HighpassShaper : public Shaper
    {
    Shaper_1Pole shaper;
    HighpassShaper( int order, const Function<Second, Frequency>& cutoff ) : shaper( order, cutoff, false ) {}
    void init( const Audio& me ) override { shaper.init( me ); }
    void reset() override { shaper.reset(); }
    Sample process(Second t, Sample s) override { return shaper.process( t, s ); }
    };

struct DelayShaper : public Shaper
    {
    Function<Second, Second> delay_length;
    Second max_delay_time = 1;
    std::vector<Sample> sample_buffer;
    Frame buffer_write_pos = 0;
    FrameRate sr;

    Sample read_interpolated(float delay_samples) const;

    DelayShaper( const Function<Second, Second>& delay_length );
    void init( const Audio& me ) override;
    void reset() override;
    Sample process(Second t, Sample s) override;
    };

struct GranularPitchShaper : public Shaper
    {
    Function<Second, float> pitch_ratio;
    Function<Second, Second> grain_length;
    Second max_grain_length;
    float max_pitch_ratio;
    std::vector<Sample> sample_buffer;
    std::array<float, 2> read_positions = {};
    std::array<Second, 2> grain_ages = {};
    Frame buffer_write_pos = 0;
    FrameRate sr;

    Sample read_interpolated(float read_position) const;

    GranularPitchShaper(
        const Function<Second, float>& pitch_ratio,
        const Function<Second, Second>& grain_length,
        Second max_grain_length = 0.1f,
        float max_pitch_ratio = 4.0f
        );
    void init(const Audio& me) override;
    void reset() override;
    Sample process(Second t, Sample s) override;
    };

}