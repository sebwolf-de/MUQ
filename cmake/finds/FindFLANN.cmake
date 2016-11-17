
find_package(PkgConfig)

if(NOT DEFINED MUQ_FLANN_DIR)
	pkg_check_modules(PC_FLANN QUIET libflann)
	set(FLANN_DEFINITIONS ${PC_FLANN_CFLAGS_OTHER})

	find_path(FLANN_INCLUDE_DIR flann/flann.hpp
          HINTS ${PC_FLANN_INCLUDEDIR} ${PC_FLANN_INCLUDE_DIRS}
          PATH_SUFFIXES flann )

	find_library(FLANN_LIBRARY NAMES flann_cpp
             HINTS ${PC_FLANN_LIBDIR} ${PC_FLANN_LIBRARY_DIRS} )

	find_library(FLANN_LIBRARY_STATIC NAMES ${library_prefix}flann_cpp_s.${static_library_suffix}
             HINTS ${PC_FLANN_LIBDIR} ${PC_FLANN_LIBRARY_DIRS} )

else()
	find_path(FLANN_INCLUDE_DIR flann/flann.hpp
	          HINTS ${MUQ_FLANN_DIR}/include
	          PATH_SUFFIXES flann NO_DEFAULT_PATH)

	find_library(FLANN_LIBRARY NAMES flann_cpp
	             HINTS ${MUQ_FLANN_DIR}/lib NO_DEFAULT_PATH)

	find_library(FLANN_LIBRARY_STATIC NAMES ${library_prefix}flann_cpp_s.${static_library_suffix}
	             HINTS ${MUQ_FLANN_DIR}/lib NO_DEFAULT_PATH)	 
endif()
set(FLANN_LIBRARIES_STATIC ${FLANN_LIBRARY_STATIC} )
set(FLANN_LIBRARIES ${FLANN_LIBRARY} )
set(FLANN_INCLUDE_DIRS ${FLANN_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(FLANN  DEFAULT_MSG
                                  FLANN_LIBRARY FLANN_INCLUDE_DIR)

mark_as_advanced(FLANN_INCLUDE_DIR FLANN_LIBRARY )
