#pragma once

#include <utility>
#include <vector>
#include <complex>

#include "flan/defines.h"
#include "flan/Utility/MP.h"

namespace flan {

class SQPVBuffer
{
public:
	struct Format 
		{ 
		Channel num_channels = 0;
		Frame num_frames = 0;
		fBin bins_per_octave = 0;
		FrameRate sample_rate = 48000; 
		std::pair<Frequency, Frequency> bandwidth;
		};

	SQPVBuffer( const SQPVBuffer & ) = delete;
	SQPVBuffer( SQPVBuffer && ) = default;
	SQPVBuffer& operator=( const SQPVBuffer & ) = delete;
	SQPVBuffer& operator=( SQPVBuffer && ) = default;
	~SQPVBuffer() = default;

	SQPVBuffer();
	SQPVBuffer( const Format & );

	MP getMP( Channel, Frame, Bin ) const;
	MP & getMP( Channel, Frame, Bin );

	const Format & get_format() const;

	fFrame time_to_frame( Second ) const;
	Second frame_to_time( fFrame ) const;
	fBin frequency_to_bin( Frequency ) const;
	Frequency bin_to_frequency( fBin ) const;
	Pitch frequencyToPitch( Frequency ) const;
	Frequency pitchToFrequency( Pitch ) const;
	UnsignedPitch binToPitch( fBin ) const;
	fBin pitchToBin( UnsignedPitch ) const;

	Channel get_num_channels() const;
	Frame get_num_frames() const;
	Bin	get_num_bins() const;
	FrameRate get_sample_rate() const;
	FrameRate get_analysis_rate() const;
	std::pair<Frequency, Frequency> getFrequencyBandwidth() const;
	std::pair<UnsignedPitch, UnsignedPitch> getPitchBandwidth() const;
	fBin getBinsPerOctave() const;
	Cycle getQ() const;
	Frequency getBinFrequency( Bin ) const;
	Second get_length() const { return frame_to_time( get_num_frames() ); }
	bool is_null() const;
	void clear_buffer();
	SQPVBuffer copy() const;
	Magnitude get_max_partial_magnitude() const;
	Magnitude get_max_partial_magnitude( Frame start_frame, Frame end_frame = 0, Bin start_bin = 0, Bin end_bin = 0 ) const;
	Frame getPeriod( Bin ) const;

	size_t get_buffer_pos( Channel, Frame, Bin ) const;

private:

	const Format format;

	// Commonly used data is computed once at construction
	const std::pair<UnsignedPitch, UnsignedPitch> pitch_bandwidth;
	const Bin num_bins;
	const Cycle Q; // "The number of cycles needed to make an analysis" - Sliding With A Constant Q
	const std::vector<Frequency> bin_frequencies;

	std::vector<MP> buffer;
};

}