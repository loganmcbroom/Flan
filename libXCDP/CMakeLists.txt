find_package( SndFile		REQUIRED ) # SndFile::sndfile
find_package( SampleRate	REQUIRED ) # SampleRate::samplerate
find_package( OpenCL		REQUIRED ) # OpenCL::OpenCL
find_package( FFTWF			REQUIRED ) # FFTWF::fftwf

add_library( libXCDP 
	src/Audio.cpp 
	src/AudioBuffer.cpp 
	src/PVOC.cpp 
	src/PVOCBuffer.cpp 
	src/Spectrum.cpp
	src/SpectrumBuffer.cpp

	src/WindowFunctions.cpp 
	src/RealFunc.cpp
	src/Utility.cpp
	src/Synthesis.cpp
	src/CLContext.cpp

	src/WDL/denormal.h
	src/WDL/heapbuf.h
	src/WDL/resample.h
	src/WDL/resample.cpp
	src/WDL/wdltypes.h

	src/CDPMapPVOC.cpp
	)

target_include_directories( libXCDP 
	PUBLIC 
		$<INSTALL_INTERFACE:include>    
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	PRIVATE
)

target_link_libraries( libXCDP 
	PRIVATE 
		SndFile::sndfile
		SampleRate::samplerate
		OpenCL::OpenCL
		FFTWF::fftwf
)	

