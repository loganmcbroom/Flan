# Find libsndfile library
# 
# Variables Set
#  	SNDFILE_FOUND 				- true if found on the system
#  	SNDFILE_LIBRARY_DIR 		- full path to library
#  	SNDFILE_INCLUDE_DIR 		- include directory path
#
# Prefer interface: SndFile::sndfile
#

find_path( SNDFILE_INCLUDE_DIR
    NAME
		sndfile.h
    PATHS
		"C:/Program Files (x86)/sndfile/include"
	PATH_SUFFIXES
		include
		libsndfile/include
  )
  
find_library( SNDFILE_LIBRARY_DIR
    NAME
		sndfile
		libsndfile-1
    PATHS
		"C:/Program Files (x86)/sndfile/lib"
	PATH_SUFFIXES
		lib
		libsndfile/lib
  )
  
include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( SndFile
	REQUIRED_VARS 
		SNDFILE_INCLUDE_DIR
		SNDFILE_LIBRARY_DIR
	HANDLE_COMPONENTS
	)

if( SNDFILE_FOUND )
	add_library( SndFile::sndfile INTERFACE IMPORTED )
	set_target_properties( SndFile::sndfile
		PROPERTIES 
			INTERFACE_INCLUDE_DIRECTORIES ${SNDFILE_INCLUDE_DIR}
			INTERFACE_LINK_LIBRARIES ${SNDFILE_LIBRARY_DIR}
	)
else()
	if( SndFile_FIND_REQUIRED )
		message( FATAL_ERROR "Couldn't find libsndfile" )
	endif()
endif()

mark_as_advanced(
    SNDFILE_INCLUDE_DIR
    SNDFILE_LIBRARY_DIR
    )
