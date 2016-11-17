
find_package(PkgConfig)

if(NOT DEFINED MUQ_HDF5_DIR)
	
	pkg_check_modules(PC_HDF5HL QUIET libhdf5_hl)
	set(HDF5HL_DEFINITIONS ${PC_HDF5HL_CFLAGS_OTHER})
	
	find_path(HDF5HL_INCLUDE_DIR hdf5_hl.h
          HINTS ${PC_HDF5HL_INCLUDEDIR} ${PC_HDF5HL_INCLUDE_DIRS}
          PATH_SUFFIXES hdf5_hl )

	find_library(HDF5HL_LIBRARY NAMES hdf5_hl
             HINTS ${PC_HDF5HL_LIBDIR} ${PC_HDF5HL_LIBRARY_DIRS} )

	find_library(HDF5HL_LIBRARY_STATIC NAMES ${library_prefix}hdf5_hl.${static_library_suffix}
             HINTS ${PC_HDF5HL_LIBDIR} ${PC_HDF5HL_LIBRARY_DIRS} )

else()
	find_path(HDF5HL_INCLUDE_DIR hdf5_hl.h
	          HINTS ${MUQ_HDF5_DIR}/include
	          PATH_SUFFIXES hdf5_hl NO_DEFAULT_PATH)

	find_library(HDF5HL_LIBRARY NAMES hdf5_hl
	             HINTS ${MUQ_HDF5_DIR}/lib NO_DEFAULT_PATH)

	find_library(HDF5HL_LIBRARY_STATIC NAMES ${library_prefix}hdf5_hl.${static_library_suffix}
	             HINTS ${MUQ_HDF5_DIR}/lib NO_DEFAULT_PATH)			 
endif()


set(HDF5HL_LIBRARIES_STATIC ${HDF5HL_LIBRARY_STATIC} )	
			 
set(HDF5HL_LIBRARIES ${HDF5HL_LIBRARY} )
set(HDF5HL_INCLUDE_DIRS ${HDF5HL_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(HDF5HL  DEFAULT_MSG
                                  HDF5HL_LIBRARY HDF5HL_INCLUDE_DIR)

mark_as_advanced(HDF5HL_INCLUDE_DIR HDF5HL_LIBRARY )