//#pragma once
//
//#include <complex>
//#include <vector>
//
//namespace flan {
//
//typedef std::complex<float> complex_f;
//
//class SpectrumBuffer
//{
//public:
//	struct Format 
//		{
//		size_t num_channels = 0, num_bins = 0;
//		size_t sample_rate = 48000;
//		};
//	
//	SpectrumBuffer();
//	SpectrumBuffer( const Format & ); 
//
//	//======================================================
//	//	Getters
//	//======================================================
//
//	complex_f getSpectra( size_t channel, size_t bin ) const;
//	
//	Format get_format() const;
//	size_t get_num_channels() const;
//	size_t get_num_bins() const;
//	size_t get_sample_rate() const;
//
//	float bin_to_frequency() const;
//	float frequency_to_bin() const;
//
//	float getMaxSpectraMagnitude() const;
//
//	//======================================================
//	//	Setters
//	//======================================================
//
//	void setSpectra( size_t channel, size_t bin, complex_f sample );
//	complex_f & getSpectra( size_t channel, size_t bin );
//	void clear_buffer();
//
//private://=================================================================================================
//
//	size_t getPos( size_t, size_t ) const;
//
//	Format format;
//	std::vector< complex_f > buffer;
//
//};
//
//}