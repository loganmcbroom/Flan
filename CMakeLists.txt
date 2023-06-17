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
	src/flan/Conversions/AudioPV.cpp
	src/flan/Conversions/AudioSPV.cpp
	src/flan/Conversions/AudioSQPV.cpp
	src/flan/Conversions/AudioGraph.cpp
	src/flan/Conversions/PVGraph.cpp

	src/flan/Audio/Audio.cpp
	src/flan/Audio/AudioBuffer.cpp 
	src/flan/Audio/AudioTexture.cpp
	src/flan/Audio/AudioConversions.cpp
	src/flan/Audio/AudioInformation.cpp
	src/flan/Audio/AudioSynthesis.cpp
	src/flan/Audio/AudioMod.cpp
	src/flan/Audio/AudioFilter.cpp
	
	src/flan/PV/PV.cpp
	src/flan/PV/PVBuffer.cpp
	src/flan/PV/PVInformation.cpp
	src/flan/PV/PVModify.cpp
	src/flan/PV/PrismFunc.h

	src/flan/SPV/SPVBuffer.cpp
	src/flan/SPV/SPV.cpp

	src/flan/SQPV/SQPVBuffer.cpp
	src/flan/SQPV/SQPV.cpp

	src/flan/Utility/Color.cpp
	src/flan/Utility/Bytes.cpp
	src/flan/Utility/Interpolator.cpp
	src/flan/Utility/View.cpp
	src/flan/Utility/Interval.cpp 
	src/flan/Utility/Rect.cpp 
	src/flan/Utility/iota_iter.h
	src/flan/Utility/buffer_access.cpp
	src/flan/Utility/execution.cpp 

	src/flan/defines.cpp 
	src/flan/WindowFunctions.cpp 
	src/flan/Function.cpp
	src/flan/FFTHelper.cpp 
	src/flan/Graph.cpp
	src/flan/Wavetable.cpp
	src/flan/DSPUtility.cpp 
	src/flan/phase_vocoder.cpp 


	# Included libs 
	src/WDL/resample.cpp

	src/r8brain/r8bbase.cpp 
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
		$<INSTALL_INTERFACE:src>    
        $<BUILD_INTERFACE: 
			${CMAKE_CURRENT_SOURCE_DIR}/src >
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

#==================================================================================================
# Installation
#==================================================================================================

include( GNUInstallDirs )
set( INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/flan )

# Copy include folder into install folder
install( DIRECTORY src DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/src )

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
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/flan
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