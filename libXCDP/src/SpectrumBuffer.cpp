#include "SpectrumBuffer.h"

#include <algorithm>

namespace xcdp {

SpectrumBuffer::SpectrumBuffer()
	: format()
	, buffer()
	{}

SpectrumBuffer::SpectrumBuffer( const Format & f )
	: format( f )
	, buffer( getNumChannels() * getNumBins() )
	{}

//======================================================
//	Getters
//======================================================

complex_d SpectrumBuffer::getSpectra( size_t channel, size_t bin ) const
	{
	return buffer[getPos( channel, bin )];
	}
SpectrumBuffer::Format SpectrumBuffer::getFormat() const
	{
	return format;
	}
size_t SpectrumBuffer::getNumChannels() const
	{
	return format.numChannels;
	}
size_t SpectrumBuffer::getNumBins() const
	{
	return format.numBins;
	}
size_t SpectrumBuffer::getSampleRate() const
	{
	return format.sampleRate;
	}
double SpectrumBuffer::binToFrequency( size_t bin ) const
	{
	return double( bin ) / double( getNumBins() ) * double( format.sampleRate );
	}
size_t SpectrumBuffer::frequencyToBin( double frequency ) const
	{
	return frequency * double( getNumBins() ) / double( format.sampleRate );
	}
double SpectrumBuffer::getMaxSpectraMagnitude() const
	{
	return std::abs(*std::max_element( buffer.begin(), buffer.end(), []( complex_d a, complex_d b )
		{
		return std::abs( a ) < std::abs( b );
		}));
	}

//======================================================
//	Setters
//======================================================

void SpectrumBuffer::setSpectra( size_t channel, size_t bin, std::complex<double> newValue )
	{
	buffer[getPos( channel, bin )] = newValue;
	}
complex_d & SpectrumBuffer::getSpectra( size_t channel, size_t bin )
	{
	return buffer[getPos( channel, bin )];
	}
void SpectrumBuffer::clearBuffer()
	{
	std::fill( buffer.begin(), buffer.end(), 0 );
	}

//======================================================
//	Private
//======================================================

size_t SpectrumBuffer::getPos( size_t channel, size_t bin ) const
	{
	return channel * getNumBins() + bin;
	}

}
