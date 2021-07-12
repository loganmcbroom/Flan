cmake_minimum_required( VERSION 3.8 )

project( flan )

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED True )
set( CMAKE_DISABLE_IN_SOURCE_BUILD True )
set( CMAKE_DEBUG_POSTFIX d )

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/" )

option( BUILD_DOC "Build documentation" ON )
if( CMAKE_BUILD_TYPE STREQUAL "Release" AND BUILD_DOC )
	find_package( Doxygen )
	if( DOXYGEN_FOUND )
		set( DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in )
		set( DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile )

		configure_file( ${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY )
		message( "Doxygen build started" ) 

		add_custom_target( doc_doxygen ALL
			COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
			COMMENT "Generating API documentation with Doxygen"
			VERBATIM )

	else( DOXYGEN_FOUND )
	  message( "Doxygen need to be installed to generate the doxygen documentation" )
	endif( DOXYGEN_FOUND )
endif()

add_library( flan 
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
	src/PVOC/PVOCModify.cpp

	src/Utility/Color.cpp
	src/Utility/Bytes.cpp
	src/Utility/Interpolator.cpp
	src/Utility/View.cpp
	src/Utility/Interval.cpp 
	src/Utility/Rect.cpp 
	#include/flan/Utility/Rect.h
	#include/flan/Utility/Interval.h
	#include/flan/Utility/vec2.h
	#include/flan/Utility/View.h

	src/WindowFunctions.cpp 
	src/Function.cpp
	src/CLContext.cpp
	src/CLProgs.cpp
	src/FFTHelper.cpp 
	src/Synthesis.cpp
	src/Graph.cpp
	src/Wavetable.cpp
	src/DSPUtility.cpp 
	
	src/WDL/resample.cpp
	include/r8brain/r8bbase.cpp 
	)
	
target_include_directories( flan 
	PUBLIC 
		$<INSTALL_INTERFACE:include>    
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	)

find_package( FFTWF REQUIRED ) # FFTWF::fftwf
target_link_libraries( flan PUBLIC FFTWF::fftwf )

option( USE_SNDFILE "Build using libsndfile" ON )
if( USE_SNDFILE )
	add_definitions( -DUSE_SNDFILE )
	find_package( SndFile REQUIRED ) # SndFile::sndfile
	target_link_libraries( flan PUBLIC SndFile::sndfile )
endif( USE_SNDFILE )
unset( USE_SNDFILE CACHE )

option( USE_OPENCL "Build using OpenCL" ON )
if( USE_OPENCL )
	add_definitions( -DUSE_OPENCL )
	find_package( OpenCL REQUIRED ) # OpenCL::OpenCL
	target_link_libraries( flan PUBLIC OpenCL::OpenCL )
endif( USE_OPENCL )
unset( USE_OPENCL CACHE )

# Copy lib to top level /lib folder after building
if( CMAKE_BUILD_TYPE STREQUAL "Release" )
    add_custom_command( TARGET flan POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:flan> ${CMAKE_CURRENT_SOURCE_DIR}/lib/flan.lib )
endif()
if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
    add_custom_command( TARGET flan POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:flan> ${CMAKE_CURRENT_SOURCE_DIR}/lib/fland.lib )
endif()

# If a testing folder is present, build it
if( EXISTS ${PROJECT_SOURCE_DIR}/tests )
	add_subdirectory( tests )
endif()