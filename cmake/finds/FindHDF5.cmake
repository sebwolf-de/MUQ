
find_package(PkgConfig)

if(NOT DEFINED MUQ_HDF5_DIR)
	pkg_check_modules(PC_HDF5 QUIET libhdf5)
	set(HDF5_DEFINITIONS ${PC_HDF5_CFLAGS_OTHER})
	
	find_path(HDF5_INCLUDE_DIR NAMES hdf5.h
          HINTS ${PC_HDF5_INCLUDEDIR} ${PC_HDF5_INCLUDE_DIRS} /usr/local/include
          PATH_SUFFIXES hdf5 )

	find_library(HDF5_LIBRARY NAMES hdf5
             HINTS ${PC_HDF5_LIBDIR} ${PC_HDF5_LIBRARY_DIRS} )

	find_library(HDF5_LIBRARY_STATIC NAMES ${library_prefix}hdf5.${static_library_suffix}
             HINTS ${PC_HDF5_LIBDIR} ${PC_HDF5_LIBRARY_DIRS} )

else()
	find_path(HDF5_INCLUDE_DIR NAMES hdf5.h
	          HINTS ${MUQ_HDF5_DIR}
		  	  PATH_SUFFIXES include NO_DEFAULT_PATH)
			  
	find_library(HDF5_LIBRARY NAMES hdf5
	             HINTS ${MUQ_HDF5_DIR}
			     PATH_SUFFIXES lib NO_DEFAULT_PATH)
	
	find_library(HDF5_LIBRARY_STATIC NAMES ${library_prefix}hdf5.${static_library_suffix}
	             HINTS ${MUQ_HDF5_DIR}
			     PATH_SUFFIXES lib NO_DEFAULT_PATH)	 
endif()

set(HDF5_LIBRARIES_STATIC ${HDF5_LIBRARY_STATIC} )	 
			 
set(HDF5_LIBRARIES ${HDF5_LIBRARY} )
set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(HDF5  DEFAULT_MSG
                                  HDF5_LIBRARY HDF5_INCLUDE_DIR)

mark_as_advanced(HDF5_INCLUDE_DIR HDF5_LIBRARY )