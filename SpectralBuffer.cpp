#include "SpectralBuffer.h"

#include <iostream>
#include <complex>

#include <fftw3.h>

#include "AudioBuffer.h"

using namespace xcdp;

SpectralBuffer::SpectralBuffer()
	: sampleRate( 48000 )
	, format( SF_FORMAT_WAV | SF_FORMAT_PCM_24 )
	{}

SpectralBuffer::SpectralBuffer( const std::string & filePath )
	: SpectralBuffer( AudioBuffer( filePath ) )
	{
	}

///This could be optimized with r2c fft
SpectralBuffer::SpectralBuffer( const AudioBuffer & in )
	{
	sampleRate = in.getSampleRate();
	format = in.info.format;

	buffer.resize( in.getNumChannels() );

	const int N = in.getNumFrames();

	std::complex<float> *inCpx = (std::complex<float>*) fftwf_malloc( sizeof(std::complex<float>) * N );
	fftwf_plan plan;

	for( int channel = 0; channel < in.getNumChannels(); ++ channel )
		{
		buffer[channel].resize( N );

		plan = fftwf_plan_dft_1d( N, (fftwf_complex*) inCpx, (fftwf_complex*) buffer[channel].data(), FFTW_FORWARD, FFTW_ESTIMATE );

		//Load a channel from in into inCpx
		for( int frame = 0; frame < N; ++frame )
			inCpx[frame] = { in.getSample( channel, frame ), 0 };

		fftwf_execute( plan );

		fftwf_destroy_plan( plan );
		}
	
	fftwf_free( inCpx ); 
	}


//======================================================
//	I/O
//======================================================

///How to handle formatting of save?
//Convert the spectrum to time-domain and save
//bool save( const std::string & filePath );

//Print some buffer data to the console
void SpectralBuffer::printSummary() const
	{
	std::cout << "\n===================================================\n";
	std::cout << "Channels: " << buffer.size() << " Bins: " << buffer[0].size() << "\n";
	std::cout << "===================================================\n\n";
	}



//======================================================
//	Getters
//======================================================

std::complex<float> SpectralBuffer::getBin( int channel, int bin ) const
	{
	return buffer[channel][bin];
	}

int SpectralBuffer::getNumChannels() const
	{
	return buffer.size();
	}

int SpectralBuffer::getNumBins() const
	{
	return buffer[0].size();
	}

float SpectralBuffer::getBinWidth() const
	{
	return float(sampleRate) / float( getNumBins() );
	}

float SpectralBuffer::getBinFreq( int bin ) const
	{
	return float( bin ) * getBinWidth();
	}

AudioBuffer SpectralBuffer::getAudio() const
	{
	//Set up output AudioBuffer, defaulting to 24bit 48000Hz WAV
	AudioBuffer out;
	out.info.format = format;
	out.info.samplerate = sampleRate;
	out.setBufferSize( getNumChannels(), getNumBins() );

	const int N = getNumBins();

	std::complex<float> * channelOut = (std::complex<float>*) fftwf_malloc( sizeof(std::complex<float>) * N );
	fftwf_plan plan;

	for( int channel = 0; channel < buffer.size(); ++ channel )
		{
		plan = fftwf_plan_dft_1d( N, (fftwf_complex*) buffer[channel].data(), (fftwf_complex*) channelOut, FFTW_BACKWARD, FFTW_ESTIMATE );

		fftwf_execute( plan );

		//Copy transformed data into out channel
		for( int frame = 0; frame < N; ++frame )
			out.setSample( channel, frame, channelOut[frame].real() / N );

		fftwf_destroy_plan( plan );
		}

	fftwf_free( channelOut );

	return out;
	}


//======================================================
//	Setters
//======================================================

void SpectralBuffer::setBin( int channel, int bin, std::complex<float> value )
	{
	buffer[channel][bin] = value;
	}


void SpectralBuffer::setNumChannels( int numChannels )
	{
	buffer.resize( numChannels );
	}


void SpectralBuffer::setNumBins( int numBins )
	{
	for( int channel = 0; channel < getNumChannels(); ++channel )
		buffer[channel].resize( numBins );
	}

void SpectralBuffer::setBufferSize( int numChannels, int numBins )
	{
	setNumChannels( numChannels );
	setNumBins( numBins );
	}

void SpectralBuffer::copyFormat( const SpectralBuffer & in )
	{
	format = in.format;
	sampleRate = in.sampleRate;
	}