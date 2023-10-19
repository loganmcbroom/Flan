//#include "flan/SpectrumBuffer.h"
//
//#include <algorithm>
//
//namespace flan {
//
//SpectrumBuffer::SpectrumBuffer()
//	: format()
//	, buffer()
//	{}
//
//SpectrumBuffer::SpectrumBuffer( const Format & f )
//	: format( f )
//	, buffer( get_num_channels() * get_num_bins() )
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
//SpectrumBuffer::Format SpectrumBuffer::get_format() const
//	{
//	return format;
//	}
//size_t SpectrumBuffer::get_num_channels() const
//	{
//	return format.num_channels;
//	}
//size_t SpectrumBuffer::get_num_bins() const
//	{
//	return format.num_bins;
//	}
//size_t SpectrumBuffer::get_sample_rate() const
//	{
//	return format.sample_rate;
//	}
//float SpectrumBuffer::bin_to_frequency() const
//	{
//	return 1.0f / frequency_to_bin();
//	}
//float SpectrumBuffer::frequency_to_bin() const
//	{
//	return float( get_num_bins() ) / float( format.sample_rate );
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
//void SpectrumBuffer::clear_buffer()
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
//	return channel * get_num_bins() + bin;
//	}
//
//}
