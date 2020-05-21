find_package( SndFile		REQUIRED ) # SndFile::sndfile
find_package( SampleRate	REQUIRED ) # SampleRate::samplerate
find_package( FFTWF			REQUIRED ) # FFTWF::fftwf

add_library( libXCDP 
	src/Audio.cpp 
	src/AudioBuffer.cpp 
	src/PVOC.cpp 
	src/PVOCBuffer.cpp 
	src/Spectrum.cpp
	src/SpectrumBuffer.cpp

	src/WindowFunctions.cpp 
	src/Function.cpp
	src/Utility.cpp
	src/Synthesis.cpp
	src/CLContext.cpp

	src/CDPMapPVOC.cpp
	src/CLProgs.cpp

	src/WDL/resample.cpp
	)

target_include_directories( libXCDP 
	PUBLIC 
		$<INSTALL_INTERFACE:include>    
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	)

target_link_libraries( libXCDP 
	PRIVATE
		SndFile::sndfile
		SampleRate::samplerate
		FFTWF::fftwf
		OpenCL::OpenCL
	)	

option( USE_OPENCL "Build using OpenCL" ON )
if( USE_OPENCL )
	add_definitions( -DUSE_OPENCL )
	find_package( OpenCL REQUIRED ) # OpenCL::OpenCL
	target_link_libraries( libXCDP PRIVATE OpenCL::OpenCL )
endif( USE_OPENCL )

