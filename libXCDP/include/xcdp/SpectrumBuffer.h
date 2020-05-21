#pragma once

#include <complex>
#include <vector>

namespace xcdp {

typedef std::complex<float> complex_f;

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

	complex_f getSpectra( size_t channel, size_t bin ) const;
	
	Format getFormat() const;
	size_t getNumChannels() const;
	size_t getNumBins() const;
	size_t getSampleRate() const;

	float binToFrequency() const;
	float frequencyToBin() const;

	float getMaxSpectraMagnitude() const;

	//======================================================
	//	Setters
	//======================================================

	void setSpectra( size_t channel, size_t bin, complex_f sample );
	complex_f & getSpectra( size_t channel, size_t bin );
	void clearBuffer();

private://=================================================================================================

	size_t getPos( size_t, size_t ) const;

	Format format;
	std::vector< complex_f > buffer;

};

}