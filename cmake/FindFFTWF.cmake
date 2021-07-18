# Find FFTW float library
#
# Variables Set
#   FFTWF_FOUND                  - true if found on the system
#   FFTWF_LIBRARY_DIR        	 - full path to library
#   FFTWF_INCLUDE_DIR            - include directory path
#
# Prefer interface: FFTWF::fftwf
#

find_path( FFTWF_INCLUDE_DIR
    NAME
		fftw3.h
	PATH_SUFFIXES
		fftw
		include
  )

find_library( FFTWF_LIBRARY_DIR
    NAME
		"libfftw3f-3"
	PATH_SUFFIXES
		fftw
		lib
)

include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( FFTWF
	REQUIRED_VARS 
		FFTWF_INCLUDE_DIR
		FFTWF_LIBRARY_DIR
	HANDLE_COMPONENTS
	)

if( FFTWF_FOUND )
	add_library( FFTWF::fftwf INTERFACE IMPORTED )
	set_target_properties( FFTWF::fftwf
		PROPERTIES 
			INTERFACE_INCLUDE_DIRECTORIES ${FFTWF_INCLUDE_DIR}
			INTERFACE_LINK_LIBRARIES ${FFTWF_LIBRARY_DIR}
	)
else()
	if( FFTWF_FIND_REQUIRED )
		message( FATAL_ERROR "Couldn't find fftwf" )
	endif()
endif()

mark_as_advanced(
        FFTWF_INCLUDE_DIR
        FFTWF_LIBRARY_DIR
        )