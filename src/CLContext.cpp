#ifdef USE_OPENCL

#include "flan/CLContext.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace flan;

CLContext::CLContext()
	{
	//std::cout << "Initializing OpenCL \n";

	//Get platform
	std::vector<cl::Platform> all_platforms;
	cl::Platform::get( &all_platforms );
	if( all_platforms.size() == 0 )
		{
		std::cout<< "No OpenCL platforms found\n";
		return;
		}
	platform = all_platforms[0];
	//std::cout << "Using platform: " << platform.getInfo<CL_PLATFORM_NAME>() << "\n";

	//Get device
	platform.getDevices( CL_DEVICE_TYPE_ALL, &devices );
	if( devices.size() == 0 )
		{
		std::cout<<" No OpenCL devices found. Check OpenCL installation!\n";
		return;
		}
	//std::cout<< "Using device: " << devices[0].getInfo<CL_DEVICE_NAME>() << "\n";
	
	context = cl::Context( { devices[0] } );

#ifdef NDEBUG
	queue = cl::CommandQueue( context, devices[0], CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE );
#else
	queue = cl::CommandQueue( context, devices[0], CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_PROFILING_ENABLE );
#endif
	
	//std::cout << "OpenCL Initialized \n\n";
	}

ProgramHelper::ProgramHelper( CLContext & context, const std::string & source )
	{
	cl_int err;

	//std::ifstream t( "OpenCL/" + file );
	//if( !t.is_open() )
	//	{
	//	std::cout << "Couldn't find " << file << std::endl;
	//	exit( -1 );
	//	}
	//std::stringstream source;
	//source << t.rdbuf();

	program = cl::Program( context.context, source, false, &err );
	if( err != CL_SUCCESS )
		{
		std::cout << "OpenCL error " << err << " while parsing:\n" << source << std::endl;
		exit( err );
		}

	err = program.build( context.devices, "-cl-denorms-are-zero" );
	if( err != CL_SUCCESS )
		{
		std::cout << "Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>( context.devices[0] ) << "\n";
		exit( err );
		}
	}


bool flan::isOpenCLAvailable()
	{
	std::vector<cl::Platform> all_platforms;
	cl::Platform::get( &all_platforms );
	return ! all_platforms.empty();
	}

#endif