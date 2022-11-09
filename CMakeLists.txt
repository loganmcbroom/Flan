cmake_minimum_required( VERSION 3.15 )

project( Flan VERSION 1.0.0 LANGUAGES CXX )

option( LOG_FUNCTION_CALLS "Log algorithm calls to std::cout" ON )
if( LOG_FUNCTION_CALLS )
	add_definitions( -Dflan_LOG_FUNCTIONS )
endif()

# Doxygen building
option( BUILD_WITH_DOCS "Build documentation" ON )
if( BUILD_WITH_DOCS )
	find_package( Doxygen )

	if( DOXYGEN_FOUND )
		message( STATUS "Doxygen found: ${DOXYGEN_EXECUTABLE} -- ${DOXYGEN_VERSION}" )
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
		message( STATUS "Doxygen needs to be installed to generate the doxygen documentation" )
	endif( DOXYGEN_FOUND )
endif()

add_library( Flan STATIC
	src/Audio/Audio.cpp
	src/Audio/AudioBuffer.cpp 
	src/Audio/AudioTexture.cpp
	src/Audio/AudioConversions.cpp
	src/Audio/AudioInformation.cpp
	src/Audio/AudioSynthesis.cpp
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
	src/Graph.cpp
	src/Wavetable.cpp
	src/DSPUtility.cpp 
	
	src/WDL/resample.cpp
	include/r8brain/r8bbase.cpp 
	)  
add_library( Flan::Flan ALIAS Flan )

target_compile_features( Flan PUBLIC cxx_std_20 )
set_target_properties( Flan PROPERTIES 
	CXX_EXTENSIONS OFF 
	DEBUG_POSTFIX d
	)

# Set includes
target_include_directories( Flan 
	PUBLIC 
		$<INSTALL_INTERFACE:include>    
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	)

#==================================================================================================
# Link Libraries
#==================================================================================================

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/" )

target_compile_definitions( Flan PUBLIC USE_SNDFILE USE_OPENCL )

find_package( FFTWF REQUIRED ) # FFTWF::fftwf
target_link_libraries( Flan PUBLIC FFTWF::fftwf )

find_package( SndFile REQUIRED ) # SndFile::sndfile
target_link_libraries( Flan PUBLIC SndFile::SndFile )

find_package( OpenCL REQUIRED ) # OpenCL::OpenCL
target_link_libraries( Flan PUBLIC OpenCL::OpenCL )

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
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/FlanConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Flan
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
export( EXPORT FlanTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/FlanTargets.cmake"
    NAMESPACE Flan::
    )

# Register package in the User Package Registry
set( CMAKE_EXPORT_PACKAGE_REGISTRY True )
export( PACKAGE Flan )

# Copy find modules into build
add_custom_command(
    TARGET Flan POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/cmake/FindFFTWF.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/cmake/ )
add_custom_command(
    TARGET Flan POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/cmake/FindSndFile.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/cmake/ )


#==================================================================================================
# Enable test folder if it exists
#==================================================================================================

if( EXISTS ${PROJECT_SOURCE_DIR}/tests )
	add_subdirectory( tests )
endif()