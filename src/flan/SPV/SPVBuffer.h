#pragma once

#include <vector>

#include "flan/defines.h"

namespace flan {

class Audio;

class SPVBuffer
{
public:
	struct Format 
		{ 
		Channel num_channels = 0;
		Frame num_frames = 0;
		Bin num_bins = 0;
		FrameRate sample_rate = 48000;
		};

	SPVBuffer( const SPVBuffer & ) = delete;
	SPVBuffer( SPVBuffer && ) = default;
	SPVBuffer& operator=( const SPVBuffer & ) = delete;
	SPVBuffer& operator=( SPVBuffer && ) = default;
	~SPVBuffer() = default;

	SPVBuffer();
	SPVBuffer( Format );

	const Format & get_format() const;
	fFrame time_to_frame( Second ) const;
	Second frame_to_time( fFrame ) const;
	fBin frequency_to_bin( Frequency ) const;
	Frequency bin_to_frequency( fBin ) const;

	Channel get_num_channels() const;
	Frame get_num_frames() const;
	Bin	get_num_bins() const;
	FrameRate get_sample_rate() const;
	FrameRate get_analysis_rate() const;

	Second get_length() const { return frame_to_time( get_num_frames() ); }
	Frequency get_height() const { return bin_to_frequency( get_num_bins() ); }
	bool is_null() const;

	MF get_MF( Channel channel, Frame frame, Bin bin ) const;
	MF & get_MF( Channel channel, Frame frame, Bin bin );
	void clear_buffer();
	SPVBuffer copy() const;

	size_t get_buffer_pos( Channel, Frame, Bin ) const;

	// Avoid using this when possible
	std::vector<MF> & get_buffer();

private:

	Format format;
	std::vector<MF> buffer;
};

}