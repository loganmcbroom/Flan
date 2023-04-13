#include "flan/Audio/Audio.h"
#include "flan/PV/PV.h"

#include <iostream>

#include "flan/FFTHelper.h"
#include "flan/WindowFunctions.h"
#include "flan/phase_vocoder.h"

using namespace flan;

PV Audio::convertToPV( Frame windowSize, Frame hopSize, Frame dftSize, flan_CANCEL_ARG_CPP ) const
	{
	flan_FUNCTION_LOG;

	const Bin numBins = dftSize / 2 + 1;
	//+1 since we analyze at start and end times
	const Frame numHops = std::ceil( getNumFrames() / hopSize ) + 1;
	
	// Set up output
	PVBuffer::Format PVFormat;
	PVFormat.numChannels = getNumChannels();
	PVFormat.numFrames = numHops;
	PVFormat.numBins = numBins;
	PVFormat.sampleRate = getSampleRate();
	PVFormat.hopSize = hopSize;
	PVFormat.windowSize = windowSize;
	PV out( PVFormat );

	// Sample Hann window
	std::vector<float> hannWindow( windowSize );
	std::transform( std::execution::par_unseq, iota_iter(0), iota_iter(windowSize), hannWindow.begin(), [&]( int i )
		{ 
		return Windows::Hann( float( i ) / float( windowSize - 1 ) ); 
		} );

	// Allocate fft buffers and phase buffer
	std::vector<float> phaseBuffer( numBins );
	FFTHelper fft( dftSize, true, false, false );

	// For each channel, do the whole thing
	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		//Set initial phase to 0
		std::fill( std::execution::par_unseq, phaseBuffer.begin(), phaseBuffer.end(), 0 );

		//For each hop, fft and save into buffer
		for( int hop = 0; hop < numHops; ++hop )
			{
			flan_CANCEL_POINT( PV() );

			// Where in the input signal we should start reading from
			const Frame inFrameStart = hopSize * hop - windowSize / 2;

			auto safeGetSample = [&]( Frame f )
				{ 
				if( f < 0 || getNumFrames() <= f ) return 0.0f;
				else return getSample( channel, f );
				};

			// Copy windowed signal into start of fft buffer
			for( Frame fftFrame = 0; fftFrame < windowSize; ++fftFrame )
				fft.getRealBuffer()[fftFrame] = safeGetSample( inFrameStart + fftFrame ) * hannWindow[fftFrame];

			// Fill the rest of the buffer with 0
			std::fill( fft.realBegin() + windowSize, fft.realEnd(), 0 );
	
			fft.r2cExecute();

			std::for_each( std::execution::par_unseq, iota_iter(0), iota_iter(numBins), [&]( Bin bin )
				{
				out.getMF( channel, hop, bin ) = phase_vocoder( phaseBuffer[bin], fft.getComplexBuffer()[bin], out.binToFrequency( bin ), frameToTime( hopSize ) );
				} );
			}
		}

	return out;
	}

PV Audio::convertToMidSidePV( Frame windowSize, Frame hop, Frame dftSize, flan_CANCEL_ARG_CPP ) const
	{
	if( getNumChannels() != 2 ) return PV();
	return convertToMidSide().convertToPV( windowSize, hop, dftSize, canceller );
	}

Audio PV::convertToAudio( flan_CANCEL_ARG_CPP ) const
	{
	flan_FUNCTION_LOG;
	
	AudioBuffer::Format audioFormat;
	audioFormat.numChannels = getNumChannels();
	audioFormat.numFrames = getNumFrames() * getHopSize();
	audioFormat.sampleRate = getSampleRate();
	Audio out( audioFormat );

	//Sample Hann window
	std::vector<float> hannWindow( getWindowSize() );
	const float windowScale = 2.67f / ( getDFTSize() * getWindowSize() / getHopSize() ); // I don't know why, converting audio->PV->audio reduces volume by an input specific amount, usually about 2.67
	std::transform( std::execution::par_unseq, iota_iter(0), iota_iter(getWindowSize()), hannWindow.begin(), [&]( int i )
		{
		return Windows::Hann( float( i ) / float( getWindowSize() - 1 ) ) * windowScale;
		} );

	std::vector<float> phaseBuffer( getNumBins() );
	FFTHelper fft( getDFTSize(), false, true, false );

	for( Channel channel = 0; channel < getNumChannels(); ++channel )
		{
		// Initial phase? Using 0, but maybe something else would work better.
		std::fill( std::execution::par_unseq, phaseBuffer.begin(), phaseBuffer.end(), 0 );

		for( Frame PVFrame = 0; PVFrame < getNumFrames(); ++PVFrame )
			{
			flan_CANCEL_POINT( Audio() );
			
			std::transform( std::execution::par_unseq, iota_iter(0), iota_iter( getNumBins() ), fft.complexBegin(), [&]( Bin bin )
				{
				return inverse_phase_vocoder( phaseBuffer[bin], getMF( channel, PVFrame, bin ), frameToTime( 1 ) );
				} );

			fft.c2rExecute();

			// Accumulate ifft output into audio buffer
			const Frame outFrameStart = getHopSize() * PVFrame - getWindowSize() / 2;
			const Frame outFrameEnd = outFrameStart + getWindowSize();
			const Frame outFrameStart_bounded = std::max( outFrameStart, 0 );
			const Frame outFrameEnd_bounded   = std::min( outFrameEnd, out.getNumFrames() );

			const Frame fftStart = outFrameStart_bounded - outFrameStart;
			const Frame fftEnd = outFrameEnd_bounded - outFrameStart;

			for( Frame fftFrame = fftStart; fftFrame < fftEnd; ++fftFrame )
				out.getSample( channel, outFrameStart + fftFrame ) += fft.getRealBuffer()[fftFrame] * hannWindow[fftFrame];
			}
		}

	return out;
	}
	
Audio PV::convertToLeftRightAudio( flan_CANCEL_ARG_CPP ) const
	{
	if( getNumChannels() != 2 ) return Audio();
	return convertToAudio( canceller ).convertToLeftRight();
	}
