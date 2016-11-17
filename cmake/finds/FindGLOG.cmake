
find_package(PkgConfig)

if(NOT DEFINED MUQ_GLOG_DIR)
	pkg_check_modules(PC_GLOG QUIET libglog)
	set(GLOG_DEFINITIONS ${PC_GLOG_CFLAGS_OTHER})

	find_path(GLOG_INCLUDE_DIR glog/logging.h
          HINTS ${PC_GLOG_INCLUDEDIR} ${PC_GLOG_INCLUDE_DIRS}
          PATH_SUFFIXES glog )

	find_library(GLOG_LIBRARY_STATIC NAMES ${library_prefix}glog.${static_library_suffix}
             HINTS ${PC_GLOG_LIBDIR} ${PC_GLOG_LIBRARY_DIRS} )
			 
	find_library(GLOG_LIBRARY NAMES glog
             HINTS ${PC_GLOG_LIBDIR} ${PC_GLOG_LIBRARY_DIRS} )

else()
	find_path(GLOG_INCLUDE_DIR glog/logging.h
	          HINTS ${MUQ_GLOG_DIR}/include
	          PATH_SUFFIXES glog NO_DEFAULT_PATH)

	find_library(GLOG_LIBRARY_STATIC NAMES ${library_prefix}glog.${static_library_suffix}
	             HINTS ${MUQ_GLOG_DIR}/lib NO_DEFAULT_PATH)
			 
	find_library(GLOG_LIBRARY NAMES glog
	             HINTS ${MUQ_GLOG_DIR}/lib NO_DEFAULT_PATH)		 
endif()

set(GLOG_LIBRARIES ${GLOG_LIBRARY} )
set(GLOG_LIBRARIES_STATIC ${GLOG_LIBRARY_STATIC} )

set(GLOG_INCLUDE_DIRS ${GLOG_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(GLOG  DEFAULT_MSG
                                  GLOG_LIBRARY GLOG_INCLUDE_DIR)

mark_as_advanced(GLOG_INCLUDE_DIR GLOG_LIBRARY )