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

	std::complex<double> *inCpx = (std::complex<double>*) fftw_malloc( sizeof(std::complex<double>) * N );
	fftw_plan plan;

	for( int channel = 0; channel < in.getNumChannels(); ++ channel )
		{
		buffer[channel].resize( N );

		plan = fftw_plan_dft_1d( N, (fftw_complex*) inCpx, (fftw_complex*) buffer[channel].data(), FFTW_FORWARD, FFTW_ESTIMATE );

		//Load a channel from in into inCpx
		for( int frame = 0; frame < N; ++frame )
			inCpx[frame] = { in.getSample( channel, frame ), 0 };

		fftw_execute( plan );

		fftw_destroy_plan( plan );
		}
	
	fftw_free( inCpx ); 
	}

std::vector<std::complex<double>> & xcdp::SpectralBuffer::operator[](int channel)
	{
	return buffer[channel];
	}

const std::vector<std::complex<double>>& xcdp::SpectralBuffer::operator[](int channel) const
	{
	return buffer[channel];
	}


//======================================================
//	I/O
//======================================================

bool xcdp::SpectralBuffer::load(const std::string& filePath)
	{
	AudioBuffer file;
	bool loaded = file.load( filePath );
	*this = SpectralBuffer( file );
	return loaded;
	}

bool SpectralBuffer::save( const std::string & filePath )
	{
	return AudioBuffer( *this ).save( filePath );
	}

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

std::complex<double> SpectralBuffer::getBin( int channel, int bin ) const
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

double SpectralBuffer::getBinWidth() const
	{
	return double(sampleRate) / double( getNumBins() );
	}

double SpectralBuffer::getBinFreq( int bin ) const
	{
	return double( bin ) * getBinWidth();
	}

AudioBuffer SpectralBuffer::getAudio() const
	{
	//Set up output AudioBuffer, defaulting to 24bit 48000Hz WAV
	AudioBuffer out;
	out.info.format = format;
	out.info.samplerate = sampleRate;
	out.setBufferSize( getNumChannels(), getNumBins() );

	const int N = getNumBins();

	std::complex<double> * channelOut = (std::complex<double>*) fftw_malloc( sizeof(std::complex<double>) * N );
	fftw_plan plan;

	for( int channel = 0; channel < buffer.size(); ++ channel )
		{
		plan = fftw_plan_dft_1d( N, (fftw_complex*) buffer[channel].data(), (fftw_complex*) channelOut, FFTW_BACKWARD, FFTW_ESTIMATE );

		fftw_execute( plan );

		//Copy transformed data into out channel
		for( int frame = 0; frame < N; ++frame )
			{
			//std::cout << channelOut[frame].real() << " " << channelOut[frame].imag() << " " << abs(channelOut[frame]) << "\n";
			out.setSample( channel, frame, channelOut[frame].real() / N );
			}

		fftw_destroy_plan( plan );
		}

	fftw_free( channelOut );

	return out;
	}


//======================================================
//	Setters
//======================================================

void SpectralBuffer::setBin( int channel, int bin, std::complex<double> value )
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

void xcdp::SpectralBuffer::setBufferSize(const SpectralBuffer& other)
	{
	setBufferSize( other.getNumChannels(), other.getNumBins() );
	}

void SpectralBuffer::copyFormat( const SpectralBuffer & in )
	{
	format = in.format;
	sampleRate = in.sampleRate;
	}