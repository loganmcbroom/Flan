#pragma once

#include <complex>
#include <mutex>

class fftwf_plan_s;

namespace flan {

// Returns the smallest power of two containing window_size
size_t power_of_2_container( size_t window_size );

struct FFTHelper
	{
	FFTHelper( uint32_t window_size, bool useR2C, bool useC2R, bool measure );
	~FFTHelper();

	void r2c_execute();
	void c2r_execute();

	size_t real_buffer_size() const { return _real_buffer_size; }
	size_t complex_buffer_size() const { return _real_buffer_size / 2 + 1; }
	float * real_begin() { return real_buffer; }
	float * real_end() { return real_buffer + real_buffer_size(); }
	std::complex<float> * complex_begin() { return complex_buffer; }
	std::complex<float> * complex_end() { return complex_buffer + complex_buffer_size(); }

	float * get_real_buffer() { return real_buffer; }
	std::complex<float> * get_complex_buffer() { return complex_buffer; }

private:
	float * real_buffer;
	std::complex<float> * complex_buffer;
	fftwf_plan_s * r2c_plan;
	fftwf_plan_s * c2r_plan;
	size_t _real_buffer_size;

	// FFTW is only thread safe for plan execution, this keeps multiple objects from creating or destroying plans at a time
	static std::recursive_mutex mutex;
	};

};