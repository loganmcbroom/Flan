#include "flan/FFTHelper.h"

#include <cassert>
#include <fftw3.h>

using namespace flan;


FFTHelper::FFTHelper( uint32_t bufferSize, bool useR2C, bool useC2R, bool measure )
	: realBuffer( fftwf_alloc_real( bufferSize ) )
	, complexBuffer( (std::complex<float>*) fftwf_alloc_complex( bufferSize / 2 + 1 ) )
	, r2cPlan( useR2C? fftwf_plan_dft_r2c_1d( bufferSize, realBuffer, (fftwf_complex*) complexBuffer, measure? FFTW_MEASURE : FFTW_ESTIMATE ) : nullptr )
	, c2rPlan( useC2R? fftwf_plan_dft_c2r_1d( bufferSize, (fftwf_complex*) complexBuffer, realBuffer, measure? FFTW_MEASURE : FFTW_ESTIMATE ) : nullptr )
	, _realBufferSize( bufferSize ) 
	{}

FFTHelper::~FFTHelper()
	{
	if( r2cPlan ) fftwf_destroy_plan( r2cPlan );
	if( c2rPlan ) fftwf_destroy_plan( c2rPlan );
	fftwf_free( realBuffer );
	fftwf_free( complexBuffer );
	}

void FFTHelper::r2cExecute() 
	{ 
	assert( r2cPlan );
	fftwf_execute( r2cPlan ); 
	}

void FFTHelper::c2rExecute() 
	{ 
	assert( c2rPlan );
	fftwf_execute( c2rPlan ); 
	}

//#else
//
//#define DJ_FFT_IMPLEMENTATION // Define this in exactly one .cpp file
//#include "flan/dj_fft.h"
//class ForwardFFTHelper
//	{
//public:
//	ForwardFFTHelper( uint32_t windowSize, float *, std::complex<float> *, fftwf_plan_s * ) 
//		: in_( windowSize )
//		, out_( windowSize, 0 )
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
//private:
//	std::vector<float> in_;
//	std::vector<std::complex<float>> out_;
//	};
//
//#endif