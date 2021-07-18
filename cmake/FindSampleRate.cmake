# Find libsamplerate library
# 
# Variables Set
#  	SAMPLERATE_FOUND 			- true if found on the system
#  	SAMPLERAT_LIBRARY_DIR 		- full path to library
#  	SAMPLERATE_INCLUDE_DIR 		- include directory path
#
# Prefer interface: SampleRate::samplerate
#

find_path( SAMPLERATE_INCLUDE_DIR
    NAME
		samplerate.h
	PATH_SUFFIXES
		libsamplerate 
		samplerate
		include
		libsamplerate/include
  )
  
find_library( SAMPLERAT_LIBRARY_DIR
    NAME
		samplerate
	PATH_SUFFIXES
		libsamplerate 
		samplerate
		lib
		libsamplerate/lib
  )
  
include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( SampleRate
	REQUIRED_VARS 
		SAMPLERATE_INCLUDE_DIR
		SAMPLERAT_LIBRARY_DIR
	HANDLE_COMPONENTS
	)

if( SAMPLERATE_FOUND )
	add_library( SampleRate::samplerate INTERFACE IMPORTED )
	set_target_properties( SampleRate::samplerate
		PROPERTIES 
			INTERFACE_INCLUDE_DIRECTORIES ${SAMPLERATE_INCLUDE_DIR}
			INTERFACE_LINK_LIBRARIES ${SAMPLERAT_LIBRARY_DIR}
	)
else()
	if( SampleRate_FIND_REQUIRED )
		message( FATAL_ERROR "Couldn't find libsamplerate" )
	endif()
endif()

mark_as_advanced(
        SAMPLERATE_INCLUDE_DIR
        SAMPLERATE_LIBRARY_DIR
        )
