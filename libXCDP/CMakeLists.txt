add_library( libXCDP 
	src/Audio/Audio.cpp
	src/Audio/AudioBuffer.cpp 
	src/Audio/AudioTexture.cpp
	src/Audio/AudioConversions.cpp
	src/Audio/AudioInformation.cpp
	src/Audio/YIN.cpp
	
	src/PVOC/PVOC.cpp
	src/PVOC/PVOCBuffer.cpp 
	src/PVOC/PVOCConversions.cpp 
	src/PVOC/PVOCInformation.cpp

	src/WindowFunctions.cpp 
	src/Function.cpp
	
	src/Synthesis.cpp
	src/Graph.cpp
	src/Wavetable.cpp
	src/DSPUtility.cpp 
	
	src/WDL/resample.cpp
	include/r8brain/r8bbase.cpp 
	src/CLContext.cpp
	src/CLProgs.cpp
	src/FFTHelper.cpp 
	  "include/xcdp/Utility/Rect.h" 
	  "include/xcdp/Utility/Interval.h" 
	  "include/xcdp/Utility/vec2.h" 
	  "src/Utility/Color.cpp" 
	  "src/Utility/Bytes.cpp" 
	  "src/Utility/Interpolator.cpp" 
	  "include/xcdp/Utility/View.h" 
	  "src/Utility/View.cpp" 
	  "src/Utility/Interval.cpp" 
	  "src/Utility/Rect.cpp" 
	  "src/PVOC/PVOCModify.cpp")

target_include_directories( libXCDP 
	PUBLIC 
		$<INSTALL_INTERFACE:include>    
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	)

find_package( FFTWF REQUIRED ) # FFTWF::fftwf
target_link_libraries( libXCDP PRIVATE FFTWF::fftwf )

option( USE_SNDFILE "Build using libsndfile" ON )
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

