# Find libsndfile library
# 
# Variables Set
#  	SndFile_FOUND 				- true if found on the system
#  	SndFile_LIBRARY_DIR 		- full path to library
#  	SndFile_INCLUDE_DIR 		- include directory path
#
# Prefer interface: SndFile::sndfile
#

find_package( PkgConfig )
pkg_check_modules( PC_SndFile QUIET SndFile )

find_path( SndFile_INCLUDE_DIR
    NAME sndfile.h
    PATHS ${PC_SndFile_INCLUDE_DIRS}
	PATH_SUFFIXES
		include
		libsndfile/include
  )
  
find_library( SndFile_LIBRARY_RELEASE
    NAME sndfile
    PATHS ${PC_SndFile_LIBRARY_DIRS}
	PATH_SUFFIXES
		lib
		libsndfile/lib
  )

find_library( SndFile_LIBRARY_DEBUG
	NAME sndfiled
	PATHS ${PC_SndFile_LIBRARY_DIRS} ${PC_SndFile_LIBRARY_DIRS}/Debug
	PATH_SUFFIXES
		lib
		libsndfile/lib
  )

# Handle Debug/Release builds
include( SelectLibraryConfigurations )
select_library_configurations( SndFile )
  
include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( SndFile
	FOUND_VAR SndFile_FOUND
	REQUIRED_VARS 
		SndFile_LIBRARY
		SndFile_INCLUDE_DIR
	VERSION_VAR SndFile_VERSION
	)

if( SndFile_FOUND )
	if( NOT TARGET SndFile::SndFile )
		add_library( SndFile::SndFile UNKNOWN IMPORTED )
	endif()

	if( SndFile_LIBRARY_RELEASE )
		set_property( TARGET SndFile::SndFile APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE )
		set_target_properties( SndFile::SndFile PROPERTIES IMPORTED_LOCATION_RELEASE "${SndFile_LIBRARY_RELEASE}" )
	endif()

	if( SndFile_LIBRARY_DEBUG )
		set_property( TARGET SndFile::SndFile APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG )
		set_target_properties( SndFile::SndFile PROPERTIES IMPORTED_LOCATION_DEBUG "${SndFile_LIBRARY_DEBUG}" )
	endif()

	set_target_properties( SndFile::SndFile PROPERTIES
		INTERFACE_COMPILE_OPTIONS "${PC_SndFile_CFLAGS_OTHER}"
		INTERFACE_INCLUDE_DIRECTORIES "${SndFile_INCLUDE_DIR}"
	)
else()
	if( QWT_FIND_REQUIRED )
		message( FATAL_ERROR "Couldn't find SndFile" )
	endif()
endif()

mark_as_advanced(
    SndFile_INCLUDE_DIR
    SndFile_LIBRARY
    )
