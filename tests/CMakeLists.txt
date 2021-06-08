add_executable( flanTest flanTest.cpp )

target_link_libraries( flanTest  
	libflan 
	winmm.lib
	)
