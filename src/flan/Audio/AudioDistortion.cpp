#include "flan/Audio/Audio.h"

#include <array>
#include <cmath>
#include <ranges>
#include "flan/Audio/filter_defs.h"

using namespace flan;
using namespace std::ranges;

Sample GainShaper::process(Second t, Sample s)
	{
	return gain(t) * s;
	}

Sample SoftClipShaper::process(Second t, Sample s)
	{
	const float drive_value = std::max( 0.0f, drive(t) );
	if( drive_value == 0.0f ) return s;
	return std::tanh( drive_value * s ) / std::tanh( drive_value );
	}

Sample RectifierShaper::process(Second, Sample s)
	{
	return full_wave ? std::abs( s ) : std::max( s, 0.0f );
	}

Sample QuantizeShaper::process(Second t, Sample s)
	{
	const Sample quantization_step = std::abs( step(t) );
	if( quantization_step == 0.0f ) return s;
	return std::round( s / quantization_step ) * quantization_step;
	}

Sample HardClipShaper::process(Second t, Sample s)
	{
	return std::clamp( gain(t) * s, -1.0f, 1.0f );
	}

Sample DCOffsetShaper::process(Second t, Sample s)
	{
	return s + offset(t);
	}

Sample SineFoldShaper::process(Second t, Sample s)
	{
	const float drive_value = std::max( 0.0f, drive(t) );
	return std::sin( 0.5f * pi * drive_value * s );
	}

Sample CubicShaper::process(Second t, Sample s)
	{
	const float drive_value = std::max( 0.0f, drive(t) );
	const float input = std::clamp( s, -1.0f, 1.0f );
	return std::clamp( input * ( 1.0f + drive_value ) - drive_value * input * input * input, -1.0f, 1.0f );
	}

Sample ExponentialShaper::process(Second t, Sample s)
	{
	const float drive_value = std::max( 0.0f, drive(t) );
	if( drive_value == 0.0f ) return s;
	const float magnitude = std::abs( s );
	const float normalization = 1.0f - std::exp( -drive_value );
	return std::copysign( ( 1.0f - std::exp( -drive_value * magnitude ) ) / normalization, s );
	}

Sample WavefolderShaper::process(Second t, Sample s)
	{
	const float input = drive(t) * s;
	const float folded = std::fmod( input + 1.0f, 4.0f ) + ( input < -1.0f ? 4.0f : 0.0f );
	return 1.0f - std::abs( folded - 2.0f );
	}

void TremoloShaper::init( const Audio& me )
	{
	sr = me.get_sample_rate();
	}

void TremoloShaper::reset()
	{
	phase = 0.0f;
	}

Sample TremoloShaper::process(Second t, Sample s)
	{
	const float mix = std::clamp( depth(t), 0.0f, 1.0f );
	const Sample modulator = 0.5f * ( 1.0f + std::sin( phase ) );
	phase = std::fmod( phase + pi2 * frequency(t) / sr, pi2 );
	if( phase < 0.0f ) phase += pi2;
	return s * ( 1.0f - mix + mix * modulator );
	}

void DownsampleShaper::init( const Audio& me )
	{
	sr = me.get_sample_rate();
	}

void DownsampleShaper::reset()
	{
	phase = 0.0f;
	held_sample = 0.0f;
	initialized = false;
	}

Sample DownsampleShaper::process(Second t, Sample s)
	{
	if( !initialized )
		{
		held_sample = s;
		initialized = true;
		}

	const FrameRate rate = std::clamp( hold_rate(t), 1.0f, sr );
	phase += rate / sr;
	if( phase >= 1.0f )
		{
		phase -= std::floor( phase );
		held_sample = s;
		}
	return held_sample;
	}

void RingModShaper::init( const Audio& me )
	{
	sr = me.get_sample_rate();
	}

void RingModShaper::reset()
	{
	phase = 0.0f;
	}

Sample RingModShaper::process(Second t, Sample s)
	{
	const float mix = std::clamp( amount(t), 0.0f, 1.0f );
	const Sample carrier = std::sin( phase );
	const Sample modulated = s * carrier;
	phase += pi2 * frequency(t) / sr;
	phase = std::fmod( phase, pi2 );
	if( phase < 0.0f ) phase += pi2;
	return s * ( 1.0f - mix ) + modulated * mix;
	}

void SlewShaper::init( const Audio& me )
	{
	sr = me.get_sample_rate();
	}

void SlewShaper::reset()
	{
	output = 0.0f;
	}

Sample SlewShaper::process(Second t, Sample s)
	{
	const Sample change = s - output;
	const Sample rate = change >= 0.0f ? std::max( 0.0f, rise_rate(t) ) : std::max( 0.0f, fall_rate(t) );
	output += std::clamp( change, -rate / sr, rate / sr );
	return output;
	}

void AllpassShaper::init( const Audio& me )
	{
	sr = me.get_sample_rate();
	const Frame size = std::max<Frame>( 2, me.time_to_frame( max_delay_time ) + 1 );
	input_buffer.resize( size );
	output_buffer.resize( size );
	reset();
	}

void AllpassShaper::reset()
	{
	std::ranges::fill( input_buffer, 0.0f );
	std::ranges::fill( output_buffer, 0.0f );
	write_pos = 0;
	}

Sample AllpassShaper::process(Second t, Sample s)
	{
	const Frame delay = std::clamp<Frame>( delay_length(t) * sr, 1, input_buffer.size() - 1 );
	const Frame read_pos = ( write_pos + input_buffer.size() - delay ) % input_buffer.size();
	const float gain = std::clamp( feedback(t), -0.999f, 0.999f );
	const Sample delayed_input = input_buffer[read_pos];
	const Sample delayed_output = output_buffer[read_pos];
	const Sample output = -gain * s + delayed_input + gain * delayed_output;
	input_buffer[write_pos] = s;
	output_buffer[write_pos] = output;
	write_pos = ( write_pos + 1 ) % input_buffer.size();
	return output;
	}

void CombShaper::init( const Audio& me )
	{
	sr = me.get_sample_rate();
	buffer.resize( std::max<Frame>( 2, me.time_to_frame( 1.0f / min_frequency ) + 1 ) );
	reset();
	}

void CombShaper::reset()
	{
	std::ranges::fill( buffer, 0.0f );
	write_pos = 0;
	}

Sample CombShaper::process(Second t, Sample s)
	{
	const Frequency frequency_value = std::max( frequency(t), min_frequency );
	const Second delay_time = 1.0f / frequency_value;
	const Frame delay = std::clamp<Frame>( delay_time * sr, 1, buffer.size() - 1 );
	const Frame read_pos = ( write_pos + buffer.size() - delay ) % buffer.size();
	const float gain = std::clamp( feedback(t), -0.999f, 0.999f );
	const Sample output = s + gain * buffer[read_pos];
	buffer[write_pos] = output;
	write_pos = ( write_pos + 1 ) % buffer.size();
	return output;
	}

void EnvelopeFollowerShaper::init( const Audio& me )
	{
	sr = me.get_sample_rate();
	}

void EnvelopeFollowerShaper::reset()
	{
	envelope = 0.0f;
	}

Sample EnvelopeFollowerShaper::process(Second t, Sample s)
	{
	const Sample target = std::abs( s );
	const Second time = target > envelope ? attack(t) : release(t);
	const float coefficient = std::exp( -1.0f / ( std::max( time, 1.0f / sr ) * sr ) );
	envelope = coefficient * envelope + ( 1.0f - coefficient ) * target;
	return s * ( 1.0f + amount(t) * envelope );
	}

void CompressorShaper::init( const Audio& me )
	{
	sr = me.get_sample_rate();
	}

void CompressorShaper::reset()
	{
	envelope = 0.0f;
	}

Sample CompressorShaper::process(Second t, Sample s)
	{
	const Sample target = std::abs( s );
	const Second time = target > envelope ? attack(t) : release(t);
	const float coefficient = std::exp( -1.0f / ( std::max( time, 1.0f / sr ) * sr ) );
	envelope = coefficient * envelope + ( 1.0f - coefficient ) * target;

	const float compression_ratio = std::max( 1.0f, ratio(t) );
	const float level = std::max( envelope, 1.0e-12f );
	const float threshold_level = std::clamp( threshold(t), 1.0e-12f, 1.0f );
	const float gain = level > threshold_level
		? std::pow( level / threshold_level, 1.0f / compression_ratio - 1.0f )
		: 1.0f;
	return s * gain;
	}

Sample DiodeShaper::process(Second t, Sample s)
	{
	const float input = s * std::max( 0.0f, drive(t) );
	const float curve = 1.0f - std::exp( -std::abs( input ) );
	const float asymmetry_value = std::clamp( asymmetry(t), -1.0f, 1.0f );
	const float positive_curve = std::clamp( curve * ( 1.0f + asymmetry_value ), 0.0f, 1.0f );
	const float negative_curve = std::clamp( curve * ( 1.0f - asymmetry_value ), 0.0f, 1.0f );
	return input >= 0.0f ? positive_curve : -negative_curve;
	}

void NoiseGateShaper::init( const Audio& me )
	{
	sr = me.get_sample_rate();
	}

void NoiseGateShaper::reset()
	{
	envelope = 0.0f;
	}

Sample NoiseGateShaper::process(Second t, Sample s)
	{
	const Sample target = std::abs( s );
	const Second time = target > envelope ? attack(t) : release(t);
	const float coefficient = std::exp( -1.0f / ( std::max( time, 1.0f / sr ) * sr ) );
	envelope = coefficient * envelope + ( 1.0f - coefficient ) * target;
	const float gate = envelope >= threshold(t) ? 1.0f : std::clamp( range(t), 0.0f, 1.0f );
	return s * gate;
	}

// Sample AsymmetricClipShaper::process(Second t, Sample s)
// 	{
// 	const float input = gain(t) * s;
// 	return input >= 0.0f
// 		? std::clamp( input, 0.0f, std::max( 0.0f, positive_limit(t) ) )
// 		: std::clamp( input, -std::max( 0.0f, negative_limit(t) ), 0.0f );
// 	}

Sample DelayShaper::read_interpolated( float delay_samples ) const 
    {
    float read_ptr = float(buffer_write_pos) - delay_samples;
    
    // Wrap read pointer into buffer bounds
    while( read_ptr < 0.0f ) read_ptr += sample_buffer.size();
    while( read_ptr >= sample_buffer.size() ) read_ptr -= sample_buffer.size();

    Frame l = Frame(read_ptr);
    Frame r = (l + 1) % sample_buffer.size();
    float frac = read_ptr - float(l);

    return sample_buffer[l] * (1.0f - frac) + sample_buffer[r] * frac;
    }

DelayShaper::DelayShaper( const Function<Second, Second>& delay_length_ )
    : delay_length( delay_length_.copy() )
    {}

void DelayShaper::init( const Audio& me ) 
    {
    sample_buffer.resize( me.time_to_frame(max_delay_time) + 2 );
    sr = me.get_sample_rate();
    }

void DelayShaper::reset()
    {
    std::ranges::fill( sample_buffer, 0.0f );
    buffer_write_pos = 0;
    }

Sample DelayShaper::process(Second t, Sample s)
    {
    sample_buffer[buffer_write_pos] = s;

    Second delay_time = delay_length(t);
    Frame delay_frames = std::clamp<Frame>(delay_time * sr, 0, sample_buffer.size() - 2 );
    Sample output = read_interpolated(delay_frames);
    buffer_write_pos = (buffer_write_pos + 1) % sample_buffer.size();

    return output;
    }

GranularPitchShaper::GranularPitchShaper(
	const Function<Second, float>& pitch_ratio_,
	const Function<Second, Second>& grain_length_,
	Second max_grain_length_,
	float max_pitch_ratio_
	)
	: pitch_ratio( pitch_ratio_.copy() )
	, grain_length( grain_length_.copy() )
	, max_grain_length( std::max( max_grain_length_, 2.0f / 48000.0f ) )
	, max_pitch_ratio( std::abs( max_pitch_ratio_ ) )
	{}

void GranularPitchShaper::init(const Audio& me)
	{
	sr = me.get_sample_rate();
	const Second max_source_time = max_grain_length * std::max( 1.0f, std::abs( max_pitch_ratio ) );
	sample_buffer.resize( std::max<Frame>( 2, me.time_to_frame( max_source_time ) + 2 ) );
	reset();
	}

void GranularPitchShaper::reset()
	{
	std::ranges::fill( sample_buffer, 0.0f );
	buffer_write_pos = 0;
	read_positions = {};
	grain_ages = { max_grain_length, max_grain_length * 0.5f };
	}

Sample GranularPitchShaper::read_interpolated(float read_position) const
	{
	while( read_position < 0.0f ) read_position += sample_buffer.size();
	while( read_position >= sample_buffer.size() ) read_position -= sample_buffer.size();

	const Frame left = Frame( read_position );
	const Frame right = ( left + 1 ) % sample_buffer.size();
	const float fraction = read_position - float( left );
	return sample_buffer[left] * ( 1.0f - fraction ) + sample_buffer[right] * fraction;
	}

Sample GranularPitchShaper::process(Second t, Sample s)
	{
	sample_buffer[buffer_write_pos] = s;

	const float ratio = std::clamp( pitch_ratio(t), -max_pitch_ratio, max_pitch_ratio );
	const Second length = std::clamp( grain_length(t), 2.0f / sr, max_grain_length );
	const float source_delay = float( length * sr * std::max( 1.0f, std::abs( ratio ) ) );
	Sample output = 0.0f;
	float window_sum = 0.0f;

	for( std::size_t grain = 0; grain < read_positions.size(); ++grain )
		{
		if( grain_ages[grain] >= length )
			{
			grain_ages[grain] = 0.0f;
			read_positions[grain] = float( buffer_write_pos ) - source_delay;
			}

		const float phase = std::clamp( float( grain_ages[grain] / length ), 0.0f, 1.0f );
		const float window = 0.5f - 0.5f * std::cos( 2.0f * float( pi ) * phase );
		output += window * read_interpolated( read_positions[grain] );
		window_sum += window;

		read_positions[grain] += ratio;
		if( read_positions[grain] < 0.0f ) read_positions[grain] += sample_buffer.size();
		if( read_positions[grain] >= sample_buffer.size() ) read_positions[grain] -= sample_buffer.size();
		grain_ages[grain] += 1.0f / sr;
		}

	buffer_write_pos = ( buffer_write_pos + 1 ) % sample_buffer.size();
	return window_sum > 0.0f ? output / window_sum : 0.0f;
	}

Audio Audio::waveshape( 
	const Function< std::pair<Second, Sample>, Sample > & shaper,
	uint16_t oversample_factor
	) const
	{
	if( is_null() ) return Audio::create_null();

	Audio oversampled = resample( get_sample_rate() * oversample_factor );
	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		runtime_execution_policy_handler( shaper.get_execution_policy(), [&]( auto policy )
			{
			std::for_each( FLAN_POLICY iota_iter( 0 ), iota_iter( oversampled.get_num_frames() ), [&]( Frame frame )
				{ 
				Sample & s = oversampled.get_sample( channel, frame );
				s = shaper( std::pair( oversampled.frame_to_time( frame ), s ) );
				} );
			} );
		}
	return oversampled.resample( get_sample_rate() );
	}

Audio Audio::waveshape_feedback( 
	const std::shared_ptr<Shaper>& shaper,
    const Function<Second, float>& feedback_amount,
	const std::vector<std::shared_ptr<Shaper>>& feedback_shapers,
    uint16_t oversample_factor
    ) const
    {
    if( is_null() ) return Audio::create_null();

    Audio oversampled = resample( get_sample_rate() * oversample_factor );

    HighpassShaper dc_blocker( 2, 15.0f );
    dc_blocker.init( oversampled );
	shaper->init( oversampled );
    for( auto& shaper : feedback_shapers )
        shaper->init( oversampled );

    for( Channel channel = 0; channel < get_num_channels(); ++channel )
        {
        dc_blocker.reset();
		shaper->reset();
        for( auto& shaper : feedback_shapers )
            shaper->reset();

        Sample feedback_sample = 0;
        
        for( Frame frame = 0; frame < oversampled.get_num_frames(); ++frame )
            { 
            Sample& s = oversampled.get_sample( channel, frame );
            const Second t = oversampled.frame_to_time( frame );

            feedback_sample = dc_blocker.process( t, feedback_sample );
            for( auto& shaper : feedback_shapers )
                feedback_sample = shaper->process( t, feedback_sample );

			s = shaper->process( t, s + feedback_amount( t ) * feedback_sample );
            feedback_sample = s;
            }
        }

    return oversampled.resample( get_sample_rate() );
    }

Audio Audio::downsample( const Function<Second, FrameRate>& new_sample_rate ) const
	{
	if( is_null() ) return Audio::create_null();

	Audio out( get_format() );
	const FrameRate source_sample_rate = get_sample_rate();
	float sample_phase = 0.0f;
	Frame held_frame = 0;

	for( Frame frame = 0; frame < get_num_frames(); ++frame )
		{
		const Second time = frame_to_time( frame );
		const FrameRate rate = std::clamp( new_sample_rate( time ), 1.0f, source_sample_rate );

		if( frame > 0 )
			{
			sample_phase += rate / source_sample_rate;
			if( sample_phase >= 1.0f )
				{
				sample_phase -= std::floor( sample_phase );
				held_frame = frame;
				}
			}

		for( Channel channel = 0; channel < get_num_channels(); ++channel )
			out.set_sample( channel, frame, get_sample( channel, held_frame ) );
		}

	return out;
	}

Audio Audio::add_moisture(
	const Function<Second, Amplitude> & amount,
	const Function<Second, Frequency> & frequency,
	const Function<Second, float> & skew,
	const Function<Second, Amplitude> & waveform
	) const
	{
	auto amount_sampled = sample_function_over_domain( amount );
	auto frequency_sampled = sample_function_over_domain( frequency );
	auto skew_sampled = sample_function_over_domain( skew );

	return waveshape( [&]( std::pair<Second, Sample> ts ) -> Sample
		{ 
		const float amount_c 	= amount_sampled[std::round(time_to_frame(ts.first))];
		const float frequency_c = frequency_sampled[std::round(time_to_frame(ts.first))];
		const float skew_c 		= skew_sampled[std::round(time_to_frame(ts.first))];

		const float power = ts.second >= 0 ? std::pow( ts.second, skew_c ) : -std::pow( -ts.second, skew_c );
		return ts.second + amount_c * ts.second * waveform( pi2 * frequency_c * power ); 
		} );
	}