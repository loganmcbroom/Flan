cmake_minimum_required( VERSION 3.8 )

project( XCDP )

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED True )

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

add_subdirectory( libXCDP )
add_subdirectory( tests )
