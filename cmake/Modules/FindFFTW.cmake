# Find FFTW double library
#
# Variables Set
#   FFTW_FOUND                  - true if found on the system
#   FFTW_LIBRARY_DIR        	- full path to library
#   FFTW_INCLUDE_DIR            - include directory path
#
# Prefer interface: FFTW::fftw
#

find_path( FFTW_INCLUDE_DIR
    NAME
		fftw3.h
	PATH_SUFFIXES
		fftw
		include
  )

find_library( FFTW_LIBRARY_DIR
    NAME
		"libfftw3-3"
	PATH_SUFFIXES
		fftw
		lib
)

include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( FFTW
	REQUIRED_VARS 
		FFTW_INCLUDE_DIR
		FFTW_LIBRARY_DIR
	HANDLE_COMPONENTS
	)

if( FFTW_FOUND )
	add_library( FFTW::fftw INTERFACE IMPORTED )
	set_target_properties( FFTW::fftw
		PROPERTIES 
			INTERFACE_INCLUDE_DIRECTORIES ${FFTW_INCLUDE_DIR}
			INTERFACE_LINK_LIBRARIES ${FFTW_LIBRARY_DIR}
	)
else()
	if( FFTW_FIND_REQUIRED )
		message( FATAL_ERROR "Couldn't find fftw" )
	endif()
endif()

mark_as_advanced(
        FFTW_INCLUDE_DIR
        FFTW_LIBRARY_DIR
        )