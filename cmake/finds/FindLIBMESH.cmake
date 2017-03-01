find_package(PkgConfig)

if (NOT DEFINED MUQ_LIBMESH_DIR)
  pkg_check_modules(PC_LIBMESH QUIET libmesh_opt libmesh_devel)

  find_path(LIBMESH_INCLUDE_DIR libmesh
    HINTS ${PC_LIBMESH_INCLUDEDIR} ${PC_LIBMESH_INCLUDE_DIRS}
    PATH_SUFFIXES libmesh )

  find_library(LIBMESH_LIBRARY NAMES mesh_opt mesh_devel
    HINTS ${PC_LIBMESH_LIBDIR} ${PC_LIBMESH_LIBRARY_DIRS} )

  find_library(LIBMESH_LIBRARY_STATIC NAMES ${library_prefix}mesh_opt.${static_library_suffix} ${library_prefix}mesh_devel.${static_library_suffix}
    HINTS ${PC_LIBMESH_LIBDIR} ${PC_LIBMESH_LIBRARY_DIRS} )

else()
  find_path(LIBMESH_INCLUDE_DIR libmesh
    HINTS ${MUQ_LIBMESH_DIR}/include
    PATH_SUFFIXES libmesh NO_DEFAULT_PATH)

  find_library(LIBMESH_OPT_LIBRARY NAMES mesh_opt
    HINTS ${MUQ_LIBMESH_DIR}/lib NO_DEFAULT_PATH)
  
  find_library(LIBMESH_DEVEL_LIBRARY NAMES mesh_devel
    HINTS ${MUQ_LIBMESH_DIR}/lib NO_DEFAULT_PATH)
  
  set(LIBMESH_LIBRARY ${LIBMESH_OPT_LIBRARY} ${LIBMESH_DEVEL_LIBRARY})

  find_library(LIBMESH_LIBRARY_STATIC NAMES ${library_prefix}mesh_opt.${static_library_suffix} ${library_prefix}mesh_devel.${static_library_suffix}
    HINTS ${MUQ_LIBMESH_DIR}/lib NO_DEFAULT_PATH)	 
endif()

set(LIBMESH_LIBRARIES_STATIC ${LIBMESH_LIBRARY_STATIC} )	
			 
set(LIBMESH_LIBRARIES ${LIBMESH_LIBRARY} )
set(LIBMESH_INCLUDE_DIRS ${LIBMESH_INCLUDE_DIR})# ${LIBMESH_MPI_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(libmesh  DEFAULT_MSG
                                  LIBMESH_LIBRARY LIBMESH_INCLUDE_DIR)

mark_as_advanced(LIBMESH_INCLUDE_DIR LIBMESH_LIBRARY )

