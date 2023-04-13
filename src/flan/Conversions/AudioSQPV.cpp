#include "flan/Audio/Audio.h"
#include "flan/SQPV/SQPV.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <deque>
#include <map>
#include <iostream>

#include "flan/phase_vocoder.h"

using namespace flan;

// Container for three complex numbers with component-wise arithmetic overloads
struct Twiddle 
	{
	std::array<std::complex<float>, 3> data;
	std::complex<float> & operator[]( size_t i ) { return data[i]; }
	const std::complex<float> & operator[]( size_t i ) const { return data[i]; }

	Twiddle operator+( const Twiddle & other ) const
		{
		return {
			data[0] + other.data[0],
			data[1] + other.data[1],
			data[2] + other.data[2],
			};
		};
	
	Twiddle operator*( const Twiddle & other ) const
		{
		return {
			data[0] * other.data[0],
			data[1] * other.data[1],
			data[2] * other.data[2],
			};
		};

	Twiddle operator+( float x ) const
		{
		return {
			data[0] + x,
			data[1] + x,
			data[2] + x,
			};
		}

	Twiddle operator-( float x ) const
		{
		return {
			data[0] - x,
			data[1] - x,
			data[2] - x,
			};
		}

	Twiddle operator*( float x ) const
		{
		return {
			data[0] * x,
			data[1] * x,
			data[2] * x,
			};
		}

	Twiddle operator/( float x ) const
		{
		return {
			data[0] / x,
			data[1] / x,
			data[2] / x,
			};
		}
	};


static inline Twiddle make_twiddle( float quality, fFrame period )
	{
	return {
		std::polar( 1.0f, pi2 * ( quality - 1 ) / period ),
		std::polar( 1.0f, pi2 * ( quality + 0 ) / period ),
		std::polar( 1.0f, pi2 * ( quality + 1 ) / period ),
		};
	}

SQPV Audio::convertToSQPV( std::pair<float, float> bandwidth, Bin binsPerOctave ) const
	{
	flan_FUNCTION_LOG;

	SQPV::Format format;
	format.numChannels = getNumChannels();
	format.numFrames = getNumFrames();
	format.binsPerOctave = binsPerOctave;
	format.bandwidth = bandwidth;
	format.sampleRate = getSampleRate();
	SQPV out( format );

	std::pair<float, float> window = std::make_pair( 0.5, -0.5 );
	const Frame longest_period = std::ceil( out.getQuality() * out.getSampleRate() / out.binToFrequency( 0 ) );
	const Frame start_offset = longest_period / 2.0f;
	
	// Compute fiddle
	const Twiddle fiddle = make_twiddle( out.getQuality(), -1.0f );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( out.getNumBins() ), [&]( Bin bin )
			{
			const size_t period = std::ceil( out.getQuality() * out.getSampleRate() / out.binToFrequency( bin ) );
			const Frame offset = -std::ceil( 1 + ( longest_period + period ) / 2.0f );
			const Twiddle twiddle = make_twiddle( out.getQuality(), period );

			const Frequency binFrequency = out.binToFrequency( bin );
			Twiddle twiddleBuffer = {0,0,0};
			Radian phaseBuffer = 0;
			for( Frame frame = 0; frame < out.getNumFrames() + start_offset; ++frame )
				{
				auto getSample_s = [&]( Frame frame ){ return frame < 0 || out.getNumFrames() <= frame ? 0.0f : getSample( channel, frame ); };
				const Sample left  = getSample_s( frame + offset + period );
				const Sample right = getSample_s( frame + offset       	  );

				twiddleBuffer = twiddle * ( twiddleBuffer + ( fiddle * left - right ) / period );

				// The sliding CQT has a delay before it starts producing output. This check waits until the output begins.
				if( frame >= start_offset )
					{
					const std::complex<float> cpx = window.first * twiddleBuffer[1] + window.second * ( twiddleBuffer[0] + twiddleBuffer[2] ) / 2.0f;
					out.getMF( channel, frame - start_offset, bin ) = phase_vocoder( phaseBuffer, cpx, binFrequency, out.frameToTime( 1 ) );
					}
				}
			} );
		}

	return out;
	}

SQPV Audio::convertToMidSideSQPV( std::pair<Frequency, Frequency> bandwidth, Bin binsPerOctave ) const
	{
	return convertToMidSide().convertToSQPV( bandwidth, binsPerOctave );
	}

Audio SQPV::convertToAudio() const
	{
	flan_FUNCTION_LOG;

	// Compute twiddles
	std::vector<Twiddle> twiddles( getNumBins() );
	for( Bin bin = 0; bin < getNumBins(); ++bin )
    	{
		const float currentPeriod = std::ceil( getQuality() * getSampleRate() / binToFrequency( bin ) );
		twiddles[bin] = make_twiddle( getQuality(), currentPeriod );
    	}

	// Setup done
	Audio::Format format;
	format.numChannels = getNumChannels();
	format.numFrames = getNumFrames();
	format.sampleRate = getSampleRate();
	Audio out( format );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		// Invert phase vocoding
		std::vector<std::complex<float>> cpxBuffer( getNumFrames() * getNumBins(), 0 );
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( getNumBins() ), [&]( Bin bin )		
			{
			Radian phaseBuffer = 0;
			for( Frame frame = 0; frame < getNumFrames(); ++frame )
				cpxBuffer[buffer_access( 0, frame, bin )] = inverse_phase_vocoder( phaseBuffer, getMF( channel, frame, bin ), frameToTime( 1 ) );
			} );

		// Invert scqt
		std::for_each( std::execution::par_unseq, iota_iter( 0 ), iota_iter( getNumFrames() ), [&]( Frame frame )
			{
			Sample sample = 0;
			for( Bin bin = 0; bin < getNumBins(); ++bin )
				sample += ( cpxBuffer[ buffer_access( 0, frame, bin ) ] * twiddles[bin][1] ).real();
			out.getSample( channel, frame ) = sample;
			} );
		}

	return out;
	}

Audio SQPV::convertToLeftRightAudio() const
	{
	return convertToAudio().convertToLeftRight();
	}