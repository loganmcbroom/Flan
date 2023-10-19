#pragma once

#include "defines.h"

namespace flan {

// Interface for time/frequency formats
class STFTInterface
{
public:
	struct Format 
		{ 
		Channel num_channels = 0;
		Frame num_frames = 0;
		Bin num_bins = 0;
		SampleRate sample_rate = 48000; 
		Frame hopSize = 0;
		Frame window_size = 0;
		//window type?
		};

	virtual fFrame time_to_frame( Second ) const = 0;
	virtual Second frame_to_time( fFrame ) const = 0;
	virtual fBin frequency_to_bin( Frequency ) const = 0;
	virtual Frequency bin_to_frequency( fBin ) const = 0;

	virtual Channel get_num_channels() const = 0;
	virtual Frame get_num_frames() const = 0;
	virtual Bin	get_num_bins( Frame f = -1 ) const = 0;
	virtual SampleRate get_sample_rate() const = 0;

	Second get_length() const { return frame_to_time( get_num_frames() ); }
	Frequency get_height() const { return bin_to_frequency( get_num_bins() ); }
	bool is_null() const 
		{ 
		return get_sample_rate() > 0
			&& get_num_channels() > 0
			&& get_num_frames() > 0
			&& get_num_bins() > 0;
		}
};

}