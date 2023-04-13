#pragma once

#include <complex>
#include <mutex>

class fftwf_plan_s;

namespace flan {

struct FFTHelper
	{
	FFTHelper( uint32_t windowSize, bool useR2C, bool useC2R, bool measure );
	~FFTHelper();

	void r2cExecute();
	void c2rExecute();

	size_t realBufferSize() const { return _realBufferSize; }
	size_t complexBufferSize() const { return _realBufferSize / 2 + 1; }
	float * realBegin() { return realBuffer; }
	float * realEnd() { return realBuffer + realBufferSize(); }
	std::complex<float> * complexBegin() { return complexBuffer; }
	std::complex<float> * complexEnd() { return complexBuffer + complexBufferSize(); }

	float * getRealBuffer() { return realBuffer; }
	std::complex<float> * getComplexBuffer() { return complexBuffer; }

private:
	float * realBuffer;
	std::complex<float> * complexBuffer;
	fftwf_plan_s * r2cPlan;
	fftwf_plan_s * c2rPlan;
	size_t _realBufferSize;

	// FFTW is only thread safe for plan execution, this keeps multiple objects from creating or destroying plans at a time
	static std::recursive_mutex mutex;
	};

};