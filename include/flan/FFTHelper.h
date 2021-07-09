#pragma once

#include <complex>

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

	float * realBuffer;
	std::complex<float> * complexBuffer;

private:
	fftwf_plan_s * r2cPlan;
	fftwf_plan_s * c2rPlan;
	size_t _realBufferSize;
	};

};