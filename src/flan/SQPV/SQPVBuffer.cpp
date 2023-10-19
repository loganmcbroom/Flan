#include "SQPVBuffer.h"

#include <algorithm>
#include <execution>

using namespace flan;

SQPVBuffer::SQPVBuffer()
	: format()
	, pitch_bandwidth()
	, num_bins()
	, Q()
	, bin_frequencies()
	, buffer()
	{}

SQPVBuffer::SQPVBuffer( const Format & f )
	: format( f )
	, pitch_bandwidth( { frequencyToPitch( getFrequencyBandwidth().first  ).p, 
						 frequencyToPitch( getFrequencyBandwidth().second ).p } )
	, num_bins( std::ceil( frequency_to_bin( getFrequencyBandwidth().second ) ) )
	, Q( 1.0f / ( std::pow( 2.0f, 1.0f / f.bins_per_octave ) - 1.0f ) )
	, bin_frequencies( [&]()
		{ 
		std::vector<Frequency> out( get_num_bins() );
		for( Bin bin = 0; bin < out.size(); ++bin )
			out[bin] = bin_to_frequency( bin ); 
		return out;
		}() )
	, buffer( get_num_channels() * get_num_frames() * get_num_bins() )
	{}

MP SQPVBuffer::getMP( Channel c, Frame f, Bin b ) const
	{
	return buffer[get_buffer_pos( c, f, b )];
	}

MP & SQPVBuffer::getMP( Channel c, Frame f, Bin b )
	{
	return buffer[get_buffer_pos( c, f, b )];
	}

const SQPVBuffer::Format & SQPVBuffer::get_format() const
	{
	return format;
	}

fFrame SQPVBuffer::time_to_frame( Second t ) const
	{
	return t * get_sample_rate();
	}

Second SQPVBuffer::frame_to_time( fFrame f ) const
	{
	return f / get_sample_rate();
	}

fBin SQPVBuffer::frequency_to_bin( Frequency f ) const
	{
	return pitchToBin( frequencyToPitch( f ).p );
	}

Frequency SQPVBuffer::bin_to_frequency( fBin bin ) const
	{
	return pitchToFrequency( Pitch( binToPitch( bin ), true ) );
	}

Pitch SQPVBuffer::frequencyToPitch( Frequency f ) const
	{
	if( f == 0 ) return std::numeric_limits<UnsignedPitch>::min();
	return Pitch( std::log2( std::abs( f ) ), f >= 0 );
	}

Frequency SQPVBuffer::pitchToFrequency( Pitch p ) const
	{
	return std::pow( 2.0f, p.p ) * ( p.positive_frequency ? 1 : -1 );
	}

UnsignedPitch SQPVBuffer::binToPitch( fBin b ) const
	{
	return b / getBinsPerOctave() + getPitchBandwidth().first;
	}

fBin SQPVBuffer::pitchToBin( UnsignedPitch p ) const
	{
	return ( p - getPitchBandwidth().first ) * getBinsPerOctave();
	}

size_t SQPVBuffer::get_buffer_pos( Channel c, Frame f , Bin b ) const
	{
	return ( c * get_num_frames() + f ) * get_num_bins() + b;
	}

Channel SQPVBuffer::get_num_channels() const
	{
	return format.num_channels;
	}

Frame SQPVBuffer::get_num_frames() const
	{
	return format.num_frames;
	}

Bin	SQPVBuffer::get_num_bins() const
	{
	return num_bins;
	}

FrameRate SQPVBuffer::get_sample_rate() const
	{
	return format.sample_rate;
	}

FrameRate SQPVBuffer::get_analysis_rate() const
	{
	return format.sample_rate;
	}

std::pair<Frequency, Frequency> SQPVBuffer::getFrequencyBandwidth() const
	{
	return format.bandwidth;
	}

std::pair<UnsignedPitch, UnsignedPitch> SQPVBuffer::getPitchBandwidth() const
	{
	return pitch_bandwidth;
	}

float SQPVBuffer::getBinsPerOctave() const
	{
	return format.bins_per_octave;
	}

Cycle SQPVBuffer::getQ() const
	{
	return Q;
	}

Frequency SQPVBuffer::getBinFrequency( Bin bin ) const
	{
	return bin_frequencies[bin];
	}

bool SQPVBuffer::is_null() const
	{
	return get_sample_rate() <= 0
		|| buffer.size() == 0
		|| getFrequencyBandwidth().first >= getFrequencyBandwidth().second;
	}

void SQPVBuffer::clear_buffer()
	{
	std::fill( std::execution::par_unseq, buffer.begin(), buffer.end(), MP{ 0, Pitch( 0, true ) } );
	}

SQPVBuffer SQPVBuffer::copy() const
	{
	SQPVBuffer out( get_format() );
	out.buffer = buffer; // Deep copy
	return out;
	}

Magnitude SQPVBuffer::get_max_partial_magnitude() const
	{
	Magnitude max_magnitude = 0;
	for( const MP & mf : buffer )
		{
		const Magnitude trueMag = std::abs( mf.m );
		if( trueMag > max_magnitude )
			max_magnitude = trueMag;
		}
	return max_magnitude;
	}

Magnitude SQPVBuffer::get_max_partial_magnitude( Frame start_frame, Frame end_frame, Bin start_bin, Bin end_bin ) const
	{
	if( end_frame == 0 ) end_frame = get_num_frames();
	if( end_bin == 0 ) end_bin = get_num_bins();

	Magnitude max_magnitude = 0;

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		for( Frame frame = start_frame; frame < end_frame; ++frame )
			for( Bin bin = start_bin; bin < end_bin; ++bin )
			    {
				const Magnitude trueMag = std::abs( getMP( channel, frame, bin ).m );
				if( trueMag > max_magnitude )
					max_magnitude = trueMag;
				}

	return max_magnitude;
	}

Frame SQPVBuffer::getPeriod( Bin bin ) const
	{
	return std::ceil( getQ() * get_sample_rate() / getBinFrequency( bin ) );
	}
