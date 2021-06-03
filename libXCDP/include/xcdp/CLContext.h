#pragma once

#ifdef USE_OPENCL

#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#include <CL/cl2.hpp>

namespace xcdp 
{

class CLContext
{
public:
	CLContext();

	cl::Platform platform;
	std::vector<cl::Device> devices;
	cl::Context context;
	cl::CommandQueue queue;
};

struct ProgramHelper
{
	ProgramHelper( CLContext & context, const std::string & source );

	cl::Program program;
};

};

#endif

