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
    const float omega = -pi2 / n;
	std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( twiddles.size() ), [&]( int i )
		{
      	twiddles[i] = std::polar( magnitude, omega * i);
		} );
	return twiddles;
	}

SPV Audio::convertToSPV( Bin numBins ) const 
	{
	flan_FUNCTION_LOG;
	
	SPV::Format format;
	format.numChannels = getNumChannels();
	format.numFrames = getNumFrames();
	format.numBins = numBins;
	format.sampleRate = getSampleRate();
	SPV out( format );

	// Alloc final output
	auto buffer_access = [&]( Frame f, Bin b ){ return f * out.getNumBins() + b; };

	// nth roots of unity for n = 2 * dftsize, clockwise order
	const auto twiddles = createTwiddles( 2 * numBins, 1.0f );
	auto getFiddle = [&]( Frame f, Bin b ) { return twiddles[ ( f * b ) % twiddles.size() ]; };

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		// Dirty buffer reuse. MF and complex<float> have the same data layout so this works and avoids a massive alloc.
		std::complex<float> * sdftBuffer = (std::complex<float> *)( out.getBuffer().data() + out.buffer_access( channel, 0, 0 ) );

		// Compute sample deltas
		std::vector<Sample> deltas( getNumFrames() );
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( deltas.size() ), [&]( Frame frame )
			{
			auto samples_s = [&]( Channel c, Frame f ){ return f >= 0 ? getSample( c, f ) : Sample( 0 ); };
			deltas[frame] = getSample( channel, frame ) - samples_s( channel, frame - 2 * out.getNumBins() );
			});

		// Compute partial sums of fiddled deltas into out (out is being reused here to avoid a buffer alloc)
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.getNumBins() ), [&]( Bin bin )
			{
			sdftBuffer[buffer_access( 0, bin )] = deltas[0];
			for( Frame frame = 1; frame < getNumFrames(); ++frame )
				sdftBuffer[buffer_access( frame, bin )] = sdftBuffer[buffer_access( frame - 1, bin )] + deltas[frame] * getFiddle( frame, bin );
			} );

		// For each frame simultaneously, a circular buffer, fiddled, is allocated with three elements
		// We then iterate through the bins of that frame computing the fiddled values into that buffer
		// Each set of three values in then convolved into the output for that frame and bin
		// The first and last bins are handled differently because there aren't three fiddled values to convolve there
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( getNumFrames() ), [&]( Frame frame )
			{
			std::complex<float> fiddled[3];
			fiddled[1] = sdftBuffer[buffer_access( frame, 0 )] * std::conj(getFiddle( frame + 1, 0 ));
			fiddled[2] = sdftBuffer[buffer_access( frame, 1 )] * std::conj(getFiddle( frame + 1, 1 ));

			// Convolve for bin = 0
			const std::complex<float> aStart = fiddled[1] + fiddled[1];
			const std::complex<float> bStart = fiddled[2].real() * 2.0f;
			const std::complex<float> convolvedStart = 0.25f * (aStart - bStart);
			sdftBuffer[buffer_access( frame, 0 )] = convolvedStart / float( out.getNumBins() * 2 );

			for( Bin bin = 1; bin < out.getNumBins() - 1; ++bin )
				{
				fiddled[(bin + 2) % 3] = sdftBuffer[buffer_access( frame, bin + 1 )] * std::conj(getFiddle( frame + 1, bin + 1 ));

				const std::complex<float> a = fiddled[(bin + 1) % 3] + fiddled[(bin + 1) % 3];
				const std::complex<float> b = fiddled[(bin + 0) % 3] + fiddled[(bin + 2) % 3];
				const std::complex<float> convolved = 0.25f * (a - b);
				sdftBuffer[buffer_access( frame, bin )] = convolved / float( out.getNumBins() * 2 );
				}

			// Convolve for bin = dftSize - 1 (the last bin)
			const std::complex<float> aEnd = fiddled[(out.getNumBins()-1 + 1) % 3] + fiddled[(out.getNumBins()-1 + 1) % 3];
			const std::complex<float> bEnd = fiddled[(out.getNumBins()-1 + 0) % 3].real() * 2.0f;
			const std::complex<float> convolvedEnd = 0.25f * (aEnd - bEnd);
			sdftBuffer[buffer_access( frame, out.getNumBins()-1 )] = convolvedEnd / float( out.getNumBins() * 2 );
			} );	
	
		// Phase vocode sdft data
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.getNumBins() ), [&]( Bin bin )		
			{
			Radian phaseBuffer = 0;
			const Frequency binFrequency = out.binToFrequency( bin );
			for( Frame frame = 0; frame < out.getNumFrames(); ++frame )
				out.getMF( channel, frame, bin ) = phase_vocoder( phaseBuffer, sdftBuffer[ frame * numBins + bin ], binFrequency, out.frameToTime( 1 ) );
			} );
		}

	return out;
	}

SPV Audio::convertToMidSideSPV( Frame dftSize ) const
	{
	return convertToMidSide().convertToSPV( dftSize );
	}

Audio SPV::convertToAudio()
	{
	flan_FUNCTION_LOG;

	Audio::Format format;
	format.numChannels = getNumChannels();
	format.numFrames = getNumFrames();
	format.sampleRate = getSampleRate();
	Audio out( format );

	std::vector<std::complex<float>> isdftIn( getNumFrames() * getNumBins(), 0 );
	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		// Invert phase vocoding
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( getNumBins() ), [&]( Bin bin )		
			{
			Radian phaseBuffer = 0;
			for( Frame frame = 0; frame < getNumFrames(); ++frame )
				isdftIn[buffer_access( 0, frame, bin)] = inverse_phase_vocoder( phaseBuffer, getMF( channel, frame, bin ), frameToTime( 1 ) );
			} );

		// Invert sdft
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.getNumFrames() ), [&]( Frame frame )
			{
			float sample = 0.0f;
			for( Bin bin = 0; bin < getNumBins(); ++bin )
				sample += isdftIn[buffer_access( 0, frame, bin)].real() * ( bin % 2 == 0 ? 1 : -1 );

			out.getSample( channel, frame ) = sample * 2.0f;
			} );
		}

	return out;
	}

Audio SPV::convertToLeftRightAudio()
	{
	return convertToAudio().convertToLeftRight();
	}
