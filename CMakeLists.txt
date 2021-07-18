cmake_minimum_required( VERSION 3.15 )

project( Flan VERSION 1.0.0 LANGUAGES CXX )

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED True )
#set( CMAKE_DISABLE_IN_SOURCE_BUILD True )
#set( CMAKE_DEBUG_POSTFIX d )

# Doxygen building
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
			VERBATIM 
			)

	else( DOXYGEN_FOUND )
	  message( "Doxygen need to be installed to generate the doxygen documentation" )
	endif( DOXYGEN_FOUND )
endif()

add_library( Flan STATIC
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
add_library( Flan::Flan ALIAS Flan )

target_compile_features( Flan PUBLIC cxx_std_17 )
set_target_properties( Flan PROPERTIES CXX_EXTENSIONS OFF )

target_include_directories( Flan 
	PUBLIC 
		$<INSTALL_INTERFACE:include>    
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	)

#==================================================================================================
# Link Libraries
#==================================================================================================

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/" )

find_package( FFTWF REQUIRED ) # FFTWF::fftwf
target_link_libraries( Flan PUBLIC FFTWF::fftwf )

#option( USE_SNDFILE "Build using libsndfile" ON )
#if( USE_SNDFILE )
	add_definitions( -DUSE_SNDFILE )
	find_package( SndFile REQUIRED ) # SndFile::sndfile
	target_link_libraries( Flan PUBLIC SndFile::sndfile )
#endif( USE_SNDFILE )
#unset( USE_SNDFILE CACHE )

#option( USE_OPENCL "Build using OpenCL" ON )
#if( USE_OPENCL )
	add_definitions( -DUSE_OPENCL )
	find_package( OpenCL REQUIRED ) # OpenCL::OpenCL
	target_link_libraries( Flan PUBLIC OpenCL::OpenCL )
#endif( USE_OPENCL )
#unset( USE_OPENCL CACHE )

#==================================================================================================
# Installation
#==================================================================================================

include( GNUInstallDirs )
set( INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/Flan )

# Copy include folder into install folder
install( DIRECTORY include/Flan/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Flan )
install( DIRECTORY include/bmp/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/bmp )

# Create export target, but don't install yet
install( 
	TARGETS Flan
    EXPORT FlanTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

# Export the export target to a *Targets.cmake script in the install folder
install( 
	EXPORT FlanTargets
    FILE FlanTargets.cmake
    NAMESPACE Flan::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Flan
)

# Create a ConfigVersion.cmake file
include( CMakePackageConfigHelpers )
write_basic_package_version_file(
    FlanConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
    )

# Generate Config file
configure_package_config_file(	
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/FlanConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/FlanConfig.cmake
    INSTALL_DESTINATION 
		${INSTALL_CONFIGDIR}
    )

# Install files: config, configversion, find modules
install( 
    FILES
		${CMAKE_CURRENT_LIST_DIR}/cmake/FindFFTWF.cmake
		${CMAKE_CURRENT_LIST_DIR}/cmake/FindSndFile.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/FlanConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/FlanConfigVersion.cmake
    DESTINATION     
        ${INSTALL_CONFIGDIR}
    )

#==================================================================================================
# Export from the build tree
#==================================================================================================

# Generate a *Targets.cmake into the build folder
export( 
	EXPORT FlanTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/FlanTargets.cmake"
    NAMESPACE Flan::
    )

# Register package in the User Package Registry
set( CMAKE_EXPORT_PACKAGE_REGISTRY True )
export( PACKAGE Flan )


#==================================================================================================
# Enable test folder if it exists
#==================================================================================================

if( EXISTS ${PROJECT_SOURCE_DIR}/tests )
	add_subdirectory( tests )
endif()