#pragma once

#include <vector>
#include <complex>

namespace xcdp {

class AudioBuffer;

class SpectralBuffer 
{
public:
	SpectralBuffer(  ); //Defaults to 24bit 48000Hz PCM WAV
	SpectralBuffer( const std::string & filePath );
	SpectralBuffer( const AudioBuffer & timeDomain );

	std::vector<std::complex<double>> & operator[]( int channel );
	const std::vector<std::complex<double>> & operator[]( int channel ) const;


	//======================================================
	//	I/O
	//======================================================

	//Load the spectrum of the given file
	bool load( const std::string & filePath );

	//Convert the spectrum to time-domain and save
	bool save( const std::string & filePath );

	//Print some buffer data to the console
	void printSummary() const;



	//======================================================
	//	Getters
	//======================================================

	//Get the sample value at the specified frame and channel
	std::complex<double> getBin( int channel, int bin ) const;

	//Get the current number of channels
	int getNumChannels() const;

	//Get the current number of frequency bins
	int getNumBins() const;

	//Get the width of each bin in Hz
	double getBinWidth() const;

	//Get the central frequency of a specific bin
	double getBinFreq( int bin ) const;

	//Convert to time domain via inverse dft
	AudioBuffer getAudio() const;


	//======================================================
	//	Setters
	//======================================================

	//Set the sample value at the specified frame and channel
	void setBin( int channel, int bin, std::complex<double> value );

	//Set the number of channels
	void setNumChannels( int numChannels );

	//Set the number of bins
	void setNumBins( int numFreqs );

	//Set both the number of channels and the number of bins, or copy the buffer size of another buffer
	void setBufferSize( int numChannels, int numBins );
	void setBufferSize( const SpectralBuffer & other );

	//Copy format and sample rate of in
	void copyFormat( const SpectralBuffer & in );



	//======================================================
	//	Members
	//======================================================

	int sampleRate; 
	int format;

	//Spectral data. First vector is channels, second is bins.
	std::vector<std::vector<std::complex<double>>> buffer;

};

} // End namespace xcdp