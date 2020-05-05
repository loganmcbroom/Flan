#pragma once

//#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#include <CL/cl2.hpp>

namespace xcdp 
{

class CLContext
{
public:
	CLContext();

	static CLContext & get();

	cl::Platform platform;
	cl::Device device;
	cl::Context context;
	cl::CommandQueue queue;
};

};

