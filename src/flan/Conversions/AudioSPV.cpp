#include "flan/Audio/Audio.h"
#include "flan/SPV/SPV.h"

#include <cmath>
#include <execution>

#include "flan/Utility/iota_iter.h"
#include "flan/defines.h"
#include "flan/phase_vocoder.h"

using namespace flan;

static std::vector<std::complex<float>> createTwiddles( size_t n, Magnitude magnitude )
	{
	std::vector<std::complex<float>> twiddles( n );
    const Radian omega = -pi2 / n;
	std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( twiddles.size() ), [&]( int i )
		{
      	twiddles[i] = std::polar( magnitude, omega * i);
		} );
	return twiddles;
	}

SPV Audio::convert_to_SPV( Bin num_bins ) const 
	{
	flan_FUNCTION_LOG;
	
	SPV::Format format;
	format.num_channels = get_num_channels();
	format.num_frames = get_num_frames();
	format.num_bins = num_bins;
	format.sample_rate = get_sample_rate();
	SPV out( format );

	auto buffer_access = [&]( Frame f, Bin b ){ return f * out.get_num_bins() + b; };

	// nth roots of unity for n = 2 * dftsize, clockwise order
	const auto twiddles = createTwiddles( 2 * num_bins, 1.0f );
	auto getFiddle = [&]( Frame f, Bin b ) { return twiddles[ ( f * b ) % twiddles.size() ]; };

	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		// Dirty buffer reuse. MF and complex<float> have the same data layout so this works and avoids a massive alloc.
		std::complex<float> * sdftBuffer = (std::complex<float> *)( out.get_buffer().data() + out.get_buffer_pos( channel, 0, 0 ) );

		// Compute sample deltas
		std::vector<Sample> deltas( get_num_frames() );
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( deltas.size() ), [&]( Frame frame )
			{
			auto samples_s = [&]( Channel c, Frame f ){ return f >= 0 ? get_sample( c, f ) : Sample( 0 ); };
			deltas[frame] = get_sample( channel, frame ) - samples_s( channel, frame - 2 * out.get_num_bins() );
			});

		// Compute partial sums of fiddled deltas into out (out is being reused here to avoid a buffer alloc)
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.get_num_bins() ), [&]( Bin bin )
			{
			sdftBuffer[buffer_access( 0, bin )] = deltas[0];
			for( Frame frame = 1; frame < get_num_frames(); ++frame )
				sdftBuffer[buffer_access( frame, bin )] = sdftBuffer[buffer_access( frame - 1, bin )] + deltas[frame] * getFiddle( frame, bin );
			} );

		// For each frame simultaneously, a circular buffer, fiddled, is allocated with three elements
		// We then iterate through the bins of that frame computing the fiddled values into that buffer
		// Each set of three values in then convolved into the output for that frame and bin
		// The first and last bins are handled differently because there aren't three fiddled values to convolve there
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( get_num_frames() ), [&]( Frame frame )
			{
			std::complex<float> fiddled[3];
			fiddled[1] = sdftBuffer[buffer_access( frame, 0 )] * std::conj(getFiddle( frame + 1, 0 ));
			fiddled[2] = sdftBuffer[buffer_access( frame, 1 )] * std::conj(getFiddle( frame + 1, 1 ));

			// Convolve for bin = 0
			const std::complex<float> aStart = fiddled[1] + fiddled[1];
			const std::complex<float> bStart = fiddled[2].real() * 2.0f;
			const std::complex<float> convolvedStart = 0.25f * (aStart - bStart);
			sdftBuffer[buffer_access( frame, 0 )] = convolvedStart / float( out.get_num_bins() * 2 );

			for( Bin bin = 1; bin < out.get_num_bins() - 1; ++bin )
				{
				fiddled[(bin + 2) % 3] = sdftBuffer[buffer_access( frame, bin + 1 )] * std::conj(getFiddle( frame + 1, bin + 1 ));

				const std::complex<float> a = fiddled[(bin + 1) % 3] + fiddled[(bin + 1) % 3];
				const std::complex<float> b = fiddled[(bin + 0) % 3] + fiddled[(bin + 2) % 3];
				const std::complex<float> convolved = 0.25f * (a - b);
				sdftBuffer[buffer_access( frame, bin )] = convolved / float( out.get_num_bins() * 2 );
				}

			// Convolve for bin = dft_size - 1 (the last bin)
			const std::complex<float> aEnd = fiddled[(out.get_num_bins()-1 + 1) % 3] + fiddled[(out.get_num_bins()-1 + 1) % 3];
			const std::complex<float> bEnd = fiddled[(out.get_num_bins()-1 + 0) % 3].real() * 2.0f;
			const std::complex<float> convolvedEnd = 0.25f * (aEnd - bEnd);
			sdftBuffer[buffer_access( frame, out.get_num_bins()-1 )] = convolvedEnd / float( out.get_num_bins() * 2 );
			} );	
	
		// Phase vocode sdft data
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.get_num_bins() ), [&]( Bin bin )		
			{
			double phase_buffer = 0;
			const Frequency binFrequency = out.bin_to_frequency( bin );
			for( Frame frame = 0; frame < out.get_num_frames(); ++frame )
				out.get_MF( channel, frame, bin ) = phase_vocoder( phase_buffer, sdftBuffer[ frame * num_bins + bin ], binFrequency, 
					out.get_analysis_rate(), out.get_sample_rate() );
			} );
		}

	return out;
	}

SPV Audio::convert_to_ms_SPV( Frame dft_size ) const
	{
	return convert_to_mid_side().convert_to_SPV( dft_size );
	}

Audio SPV::convert_to_audio()
	{
	flan_FUNCTION_LOG;

	Audio::Format format;
	format.num_channels = get_num_channels();
	format.num_frames = get_num_frames();
	format.sample_rate = get_sample_rate();
	Audio out( format );

	std::vector<std::complex<float>> isdftIn( get_num_frames() * get_num_bins(), 0 );
	for( Channel channel = 0; channel < get_num_channels(); ++channel )
		{
		// Invert phase vocoding
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( get_num_bins() ), [&]( Bin bin )		
			{
			double phase_buffer = 0;
			for( Frame frame = 0; frame < get_num_frames(); ++frame )
				isdftIn[get_buffer_pos( 0, frame, bin)] = inverse_phase_vocoder( phase_buffer, get_MF( channel, frame, bin ), get_analysis_rate() );
			} );

		// Invert sdft
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.get_num_frames() ), [&]( Frame frame )
			{
			float sample = 0.0f;
			for( Bin bin = 0; bin < get_num_bins(); ++bin )
				sample += isdftIn[get_buffer_pos( 0, frame, bin)].real() * ( bin % 2 == 0 ? 1 : -1 );

			out.get_sample( channel, frame ) = sample * 2.0f;
			} );
		}

	return out;
	}

Audio SPV::convert_to_lr_audio()
	{
	return convert_to_audio().convert_to_left_right();
	}
