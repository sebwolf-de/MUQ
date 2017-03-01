include(ExternalProject)

#NOTE: we need to have already called BuildHDF5 for this to work
ExternalProject_Get_Property( HDF5 source_dir )
ExternalProject_Get_Property( HDF5 binary_dir )
		   
set( HDF5HL_INCLUDE_DIRS "${HDF5_INSTALL_DIR}include" )
if(MUQ_USE_OPENMPI)
	set( HDF5HL_LIBRARIES ${HDF5_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl${CMAKE_STATIC_LIBRARY_SUFFIX})
else()
	set( HDF5HL_LIBRARIES ${HDF5_INSTALL_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}hdf5_hl${CMAKE_SHARED_LIBRARY_SUFFIX})
endif()
set( HDF5HL_LIBRARIES_STATIC ${HDF5_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl${CMAKE_STATIC_LIBRARY_SUFFIX})
