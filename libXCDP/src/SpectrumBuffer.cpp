//#include "xcdp/SpectrumBuffer.h"
//
//#include <algorithm>
//
//namespace xcdp {
//
//SpectrumBuffer::SpectrumBuffer()
//	: format()
//	, buffer()
//	{}
//
//SpectrumBuffer::SpectrumBuffer( const Format & f )
//	: format( f )
//	, buffer( getNumChannels() * getNumBins() )
//	{}
//
////======================================================
////	Getters
////======================================================
//
//complex_f SpectrumBuffer::getSpectra( size_t channel, size_t bin ) const
//	{
//	return buffer[getPos( channel, bin )];
//	}
//SpectrumBuffer::Format SpectrumBuffer::getFormat() const
//	{
//	return format;
//	}
//size_t SpectrumBuffer::getNumChannels() const
//	{
//	return format.numChannels;
//	}
//size_t SpectrumBuffer::getNumBins() const
//	{
//	return format.numBins;
//	}
//size_t SpectrumBuffer::getSampleRate() const
//	{
//	return format.sampleRate;
//	}
//float SpectrumBuffer::binToFrequency() const
//	{
//	return 1.0f / frequencyToBin();
//	}
//float SpectrumBuffer::frequencyToBin() const
//	{
//	return float( getNumBins() ) / float( format.sampleRate );
//	}
//float SpectrumBuffer::getMaxSpectraMagnitude() const
//	{
//	return std::abs(*std::max_element( buffer.begin(), buffer.end(), []( complex_f a, complex_f b )
//		{
//		return std::abs( a ) < std::abs( b );
//		}));
//	}
//
////======================================================
////	Setters
////======================================================
//
//void SpectrumBuffer::setSpectra( size_t channel, size_t bin, complex_f newValue )
//	{
//	buffer[getPos( channel, bin )] = newValue;
//	}
//complex_f & SpectrumBuffer::getSpectra( size_t channel, size_t bin )
//	{
//	return buffer[getPos( channel, bin )];
//	}
//void SpectrumBuffer::clearBuffer()
//	{
//	std::fill( buffer.begin(), buffer.end(), 0 );
//	}
//
////======================================================
////	Private
////======================================================
//
//size_t SpectrumBuffer::getPos( size_t channel, size_t bin ) const
//	{
//	return channel * getNumBins() + bin;
//	}
//
//}
