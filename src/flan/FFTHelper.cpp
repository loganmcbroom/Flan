#include "flan/FFTHelper.h"

#include <cassert>
#include <fftw3.h>
#include "FFTHelper.h"

using namespace flan;

std::recursive_mutex FFTHelper::mutex;

size_t flan::power_of_2_container( size_t window_size )
	{
	return std::pow( 2, (int) std::ceil( std::log2( window_size ) ) );
	}

FFTHelper::FFTHelper( uint32_t buffer_size, bool useR2C, bool useC2R, bool measure )
	{
	// We need to lock a mutex here so we can't use initializer syntax
	std::lock_guard<std::recursive_mutex> lock( mutex );

	real_buffer = fftwf_alloc_real( buffer_size );
	complex_buffer = (std::complex<float>*) fftwf_alloc_complex( buffer_size / 2 + 1 );
	r2c_plan = useR2C? fftwf_plan_dft_r2c_1d( buffer_size, real_buffer, (fftwf_complex*) complex_buffer, measure? FFTW_MEASURE : FFTW_ESTIMATE ) : nullptr;
	c2r_plan = useC2R? fftwf_plan_dft_c2r_1d( buffer_size, (fftwf_complex*) complex_buffer, real_buffer, measure? FFTW_MEASURE : FFTW_ESTIMATE ) : nullptr;
	_real_buffer_size = buffer_size;
	}

FFTHelper::~FFTHelper()
	{
	std::lock_guard<std::recursive_mutex> lock( mutex );

	if( r2c_plan ) fftwf_destroy_plan( r2c_plan );
	if( c2r_plan ) fftwf_destroy_plan( c2r_plan );
	fftwf_free( real_buffer );
	fftwf_free( complex_buffer );
	}

void FFTHelper::r2c_execute() 
	{ 
	assert( r2c_plan );
	fftwf_execute( r2c_plan ); 
	}

void FFTHelper::c2r_execute() 
	{ 
	assert( c2r_plan );
	fftwf_execute( c2r_plan ); 
	}

	// #else
	//
	// #define DJ_FFT_IMPLEMENTATION // Define this in exactly one .cpp file
	// #include "flan/dj_fft.h"
	// class ForwardFFTHelper
	//	{
	// public:
	//	ForwardFFTHelper( uint32_t window_size, float *, std::complex<float> *, fftwf_plan_s * )
	//		: in_( window_size )
	//		, out_( window_size, 0 )
	//		{}
	//
	//	void execute()
	//		{
	//		std::vector<std::complex<float>> in_c( in_.size() );
	//		for( int i = 0; i < in_.size(); ++i ) in_c[i] = in_[i]; //convert to complex
	//		out_ = dj::fft1d_gpu( in_c, dj::fft_dir::DIR_FWD );
	//		}
	//
	//	std::complex<float>  * outBuffer() { return out_.data(); }
	//
	//	float in( uint32_t n ) { return in_[n]; }
	//	void setIn( uint32_t n, float v ) { in_[n] = v; }
	//
	//	std::complex<float> out( uint32_t n ) { return out_[n]; }
	//	void setOut( uint32_t n, std::complex<float> v ) { out_[n] = v; }
	//
	// private:
	//	std::vector<float> in_;
	//	std::vector<std::complex<float>> out_;
	//	};
	//
	// #endif
