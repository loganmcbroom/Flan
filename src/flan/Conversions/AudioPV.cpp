#include "flan/Audio/Audio.h"
#include "flan/PV/PV.h"

#include <iostream>

#include "flan/FFTHelper.h"
#include "flan/WindowFunctions.h"
#include "flan/phase_vocoder.h"

using namespace flan;

PV Audio::convert_to_PV( Frame window_size, Frame hopSize, Frame dft_size, flan_CANCEL_ARG_CPP ) const
	{
	flan_FUNCTION_LOG;

	const Bin num_bins = dft_size / 2 + 1;
	//+1 since we analyze at start and end times
	const Frame numHops = std::ceil( get_num_frames() / hopSize ) + 1;
	
	// Set up output
	PVBuffer::Format PVFormat;
	PVFormat.num_channels = get_num_channels();
	PVFormat.num_frames = numHops;
	PVFormat.num_bins = num_bins;
	PVFormat.sample_rate = get_sample_rate();
	PVFormat.analysis_rate = get_sample_rate() / hopSize;
	PVFormat.window_size = window_size;
	PV out( PVFormat );

	// Sample hann window
	std::vector<float> hann_window( window_size );
	std::transform( std::execution::par_unseq, iota_iter(0), iota_iter(window_size), hann_window.begin(), [&]( int i )
		{ 
		return Windows::hann( float( i ) / float( window_size - 1 ) ); 
		} );

	// Allocate fft buffers and phase buffer
	std::vector<double> phase_buffer( num_bins );
	FFTHelper fft( dft_size, true, false, false );

	// For each channel, do the whole thing
	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		//Set initial phase to 0
		std::fill( std::execution::par_unseq, phase_buffer.begin(), phase_buffer.end(), 0 );

		//For each hop, fft and save into buffer
		for( Frame pvFrame = 0; pvFrame < numHops; ++pvFrame )
			{
			flan_CANCEL_POINT( PV() );

			// Where in the input signal we should start reading from
			const Frame in_frameStart = hopSize * pvFrame - window_size / 2;

			auto safe_get_sample = [&]( Frame f )
				{ 
				if( f < 0 || get_num_frames() <= f ) return 0.0f;
				else return get_sample( channel, f );
				};

			// Copy windowed signal into start of fft buffer
			for( Frame fftFrame = 0; fftFrame < window_size; ++fftFrame )
				fft.get_real_buffer()[fftFrame] = safe_get_sample( in_frameStart + fftFrame ) * hann_window[fftFrame];

			// Fill the rest of the buffer with 0
			std::fill( fft.real_begin() + window_size, fft.real_end(), 0 );
	
			fft.r2c_execute();

			std::for_each( std::execution::par_unseq, iota_iter(0), iota_iter(num_bins), [&]( Bin bin )
				{
				out.get_MF( channel, pvFrame, bin ) = phase_vocoder( phase_buffer[bin], fft.get_complex_buffer()[bin], 
					out.bin_to_frequency( bin ), out.get_analysis_rate(), out.get_sample_rate() );
				} );
			}
		}

	return out;
	}

PV Audio::convert_to_ms_PV( Frame window_size, Frame hop, Frame dft_size, flan_CANCEL_ARG_CPP ) const
	{
	if( get_num_channels() != 2 ) return PV();
	return convert_to_mid_side().convert_to_PV( window_size, hop, dft_size, canceller );
	}

Audio PV::convert_to_audio( flan_CANCEL_ARG_CPP ) const
	{
	flan_FUNCTION_LOG;
	
	AudioBuffer::Format audio_format;
	audio_format.num_channels = get_num_channels();
	audio_format.num_frames = get_num_frames() * get_hop_size();
	audio_format.sample_rate = get_sample_rate();
	Audio out( audio_format );

	//Sample hann window
	std::vector<float> hann_window( get_window_size() );
	const float window_scale = 2.67f / ( get_dft_size() * get_window_size() / get_hop_size() ); // I don't know why, converting audio->PV->audio reduces volume by an input specific amount, usually about 2.67
	std::transform( std::execution::par_unseq, iota_iter(0), iota_iter(get_window_size()), hann_window.begin(), [&]( int i )
		{
		return Windows::hann( float( i ) / float( get_window_size() - 1 ) ) * window_scale;
		} );

	std::vector<double> phase_buffer( get_num_bins() );
	FFTHelper fft( get_dft_size(), false, true, false );

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		// Initial phase? Using 0, but maybe something else would work better.
		std::fill( std::execution::par_unseq, phase_buffer.begin(), phase_buffer.end(), 0 );

		for( Frame pv_frame = 0; pv_frame < get_num_frames(); ++pv_frame )
			{
			flan_CANCEL_POINT( Audio() );
			
			std::transform( std::execution::par_unseq, iota_iter(0), iota_iter( get_num_bins() ), fft.complex_begin(), [&]( Bin bin )
				{
				return inverse_phase_vocoder( phase_buffer[bin], get_MF( channel, pv_frame, bin ), get_analysis_rate() );
				} );

			fft.c2r_execute();

			// Accumulate ifft output into audio buffer
			const Frame out_frameStart = get_hop_size() * pv_frame - get_window_size() / 2;
			const Frame out_frameEnd = out_frameStart + get_window_size();
			const Frame out_frameStart_bounded = std::max( out_frameStart, 0 );
			const Frame out_frameEnd_bounded   = std::min( out_frameEnd, out.get_num_frames() );

			const Frame fftStart = out_frameStart_bounded - out_frameStart;
			const Frame fftEnd = out_frameEnd_bounded - out_frameStart;

			for( Frame fftFrame = fftStart; fftFrame < fftEnd; ++fftFrame )
				out.get_sample( channel, out_frameStart + fftFrame ) += fft.get_real_buffer()[fftFrame] * hann_window[fftFrame];
			}
		}

	return out;
	}
	
Audio PV::convert_to_lr_audio( flan_CANCEL_ARG_CPP ) const
	{
	if( get_num_channels() != 2 ) return Audio();
	return convert_to_audio( canceller ).convert_to_left_right();
	}
