#include "SPVBuffer.h"

#include "flan/Utility/execution.h"

using namespace flan;

SPVBuffer::SPVBuffer()
	: format()
	, buffer()
	{}

SPVBuffer::SPVBuffer( Format _format )
	: format( _format )
	, buffer( get_num_channels() * get_num_frames() * get_num_bins() )
	{}

const SPVBuffer::Format & SPVBuffer::get_format() const 
	{
	return format;
	}

fFrame SPVBuffer::time_to_frame( Second t ) const 
	{
	return t * get_sample_rate();
	}

Second SPVBuffer::frame_to_time( fFrame f ) const 
	{
	return f / get_sample_rate();
	}

fBin SPVBuffer::frequency_to_bin( Frequency f ) const 
	{
	return f * get_num_bins() / get_sample_rate();
	}

Frequency SPVBuffer::bin_to_frequency( fBin b ) const 
	{
	return b * get_sample_rate() / get_num_bins();
	}

Channel SPVBuffer::get_num_channels() const
	{
	return format.num_channels;
	}

Frame SPVBuffer::get_num_frames() const
	{
	return format.num_frames;
	}

Bin	SPVBuffer::get_num_bins() const
	{
	return format.num_bins;
	}

FrameRate SPVBuffer::get_sample_rate() const
	{
	return format.sample_rate;
	}

FrameRate SPVBuffer::get_analysis_rate() const
	{
	return format.sample_rate;
	}

bool SPVBuffer::is_null() const 
	{ 
	return get_sample_rate() <= 0
		|| buffer.size() == 0;
	}

MF SPVBuffer::get_MF( Channel channel, Frame frame, Bin bin ) const 
	{
	return buffer[ get_buffer_pos( channel, frame, bin ) ];
	}

MF & SPVBuffer::get_MF( Channel channel, Frame frame, Bin bin )
	{
	return buffer[ get_buffer_pos( channel, frame, bin ) ];
	}

void SPVBuffer::clear_buffer()
	{
	std::fill( FLAN_PAR_UNSEQ buffer.begin(), buffer.end(), MF{0,0} );
	}
	
SPVBuffer SPVBuffer::copy() const
	{
	SPVBuffer out;
	out.format = format;
	out.buffer = buffer; // Deep copy
	return out;
	}

size_t SPVBuffer::get_buffer_pos( Channel c, Frame f, Bin b ) const
	{
	return ( c * get_num_frames() + f ) * get_num_bins() + b;
	}

std::vector<MF> & SPVBuffer::get_buffer()
	{
	return buffer;
	}


