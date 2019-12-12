#pragma once

#include <complex>
#include <vector>

namespace xcdp {

typedef std::complex<double> complex_d;

class SpectrumBuffer
{
public:
	struct Format 
		{
		size_t numChannels = 0, numBins = 0;
		size_t sampleRate = 48000;
		};
	
	SpectrumBuffer();
	SpectrumBuffer( const Format & ); 

	//======================================================
	//	Getters
	//======================================================

	complex_d getSpectra( size_t channel, size_t bin ) const;
	
	Format getFormat() const;
	size_t getNumChannels() const;
	size_t getNumBins() const;

	double binToFrequency( size_t sample ) const;
	size_t frequencyToBin( double time ) const;

	double getMaxSpectraMagnitude() const;

	//======================================================
	//	Setters
	//======================================================

	void setSpectra( size_t channel, size_t bin, std::complex<double> sample );
	complex_d & getSpectra( size_t channel, size_t bin );
	void clearBuffer();

private://=================================================================================================

	size_t getPos( size_t, size_t ) const;

	Format format;
	std::vector< complex_d > buffer;

};

}