add_executable( XCDPTest XCDPTest.cpp )

target_link_libraries( XCDPTest  
	libXCDP 
	winmm.lib
	)

add_definitions( -DUSE_OPENCL )



