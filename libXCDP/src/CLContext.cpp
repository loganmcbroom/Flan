#include "CLContext.h"

#include <iostream>

using namespace xcdp;

CLContext::CLContext()
	{
	std::cout << "Initializing OpenCL \n";

	//Get platform
	std::vector<cl::Platform> all_platforms;
	cl::Platform::get( &all_platforms );
	if( all_platforms.size() == 0 )
		{
		std::cout<< "No platforms found\n";
		std::exit( 1 );
		}
	platform = all_platforms[0];
	std::cout << "Using platform: " << platform.getInfo<CL_PLATFORM_NAME>() << "\n";

	//Get device
	std::vector<cl::Device> all_devices;
	platform.getDevices( CL_DEVICE_TYPE_ALL, &all_devices );
	if( all_devices.size() == 0 )
		{
		std::cout<<" No devices found. Check OpenCL installation!\n";
		exit(1);
		}
	device = all_devices[0];
	std::cout<< "Using device: " << device.getInfo<CL_DEVICE_NAME>() << "\n";
	
	context = cl::Context( { device } );

#ifdef NDEBUG
	queue = cl::CommandQueue( context, device, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE );
#else
	queue = cl::CommandQueue( context, device, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_PROFILING_ENABLE );
#endif
	
	std::cout << "OpenCL Initialized \n\n";
	}

CLContext & CLContext::get()
	{
	static CLContext instance;
	return instance;
	}