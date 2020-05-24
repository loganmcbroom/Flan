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

option( USE_FFTW "Build using FFTW" ON )
if( USE_FFTW )
	add_definitions( -DUSE_FFTW )
	find_package( FFTWF REQUIRED ) # FFTWF::fftwf
	target_link_libraries( libXCDP PRIVATE FFTWF::fftwf )
endif( USE_FFTW )
unset( USE_FFTW CACHE )

option( USE_SNDFILE "Build using libsndfile" Off )
if( USE_SNDFILE )
	add_definitions( -DUSE_SNDFILE )
	find_package( SndFile REQUIRED ) # SndFile::sndfile
	target_link_libraries( libXCDP PRIVATE SndFile::sndfile )
endif( USE_SNDFILE )
unset( USE_SNDFILE CACHE )

option( USE_OPENCL "Build using OpenCL" ON )
if( USE_OPENCL )
	add_definitions( -DUSE_OPENCL )
	find_package( OpenCL REQUIRED ) # OpenCL::OpenCL
	target_link_libraries( libXCDP PRIVATE OpenCL::OpenCL )
endif( USE_OPENCL )
unset( USE_OPENCL CACHE )

