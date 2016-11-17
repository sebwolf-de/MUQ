FIND_PACKAGE(Trilinos)

if(MUQ_SACADO_DIR AND NOT DEFINED MUQ_TRILINOS_DIR)
  set(MUQ_TRILINOS_DIR ${MUQ_SACADO_DIR})
endif()

if(NOT DEFINED MUQ_TRILINOS_DIR)    
    find_path(SACADO_INCLUDE_DIR Sacado.hpp
              HINTS ${Trilinos_DIR} ${Trilinos_INCLUDE_DIRS})
	
    find_library(SACADO_LIBRARY NAMES sacado
                 HINTS ${Trilinos_DIR} ${Trilinos_LIBRARY_DIRS})
		
    find_library(SACADO_LIBRARY_STATIC NAMES ${library_prefix}sacado.${static_library_suffix}
                 HINTS ${Trilinos_DIR} ${Trilinos_LIBRARY_DIRS})

    find_path(TEUCHOS_INCLUDE_DIR Teuchos_config.h
              HINTS ${MUQ_TRILINOS_DIR}/include NO_DEFAULT_PATH)

    find_library(TEUCHOS_LIBRARY NAMES teuchoscore
                 HINTS ${MUQ_TRILINOS_DIR}/lib NO_DEFAULT_PATH)

    find_library(TEUCHOS_LIBRARY_STATIC NAMES ${library_prefix}teuchoscore.${static_library_suffix}
                 HINTS ${MUQ_TRILINOS_DIR}/lib NO_DEFAULT_PATH)

else()
    
    find_path(SACADO_INCLUDE_DIR Sacado.hpp
              HINTS ${MUQ_TRILINOS_DIR}/include NO_DEFAULT_PATH)
 
    find_library(SACADO_LIBRARY NAMES sacado
                 HINTS ${MUQ_TRILINOS_DIR}/lib NO_DEFAULT_PATH)
		
    find_library(SACADO_LIBRARY_STATIC NAMES ${library_prefix}sacado.${static_library_suffix}
                 HINTS ${MUQ_TRILINOS_DIR}/lib NO_DEFAULT_PATH)

    find_path(TEUCHOS_INCLUDE_DIR Teuchos_config.h
              HINTS ${MUQ_TRILINOS_DIR}/include NO_DEFAULT_PATH)

    find_library(TEUCHOS_LIBRARY NAMES teuchoscore
                 HINTS ${MUQ_TRILINOS_DIR}/lib NO_DEFAULT_PATH)

    find_library(TEUCHOS_LIBRARY_STATIC NAMES ${library_prefix}teuchoscore.${static_library_suffix}
                 HINTS ${MUQ_TRILINOS_DIR}/lib NO_DEFAULT_PATH)

endif()

if(TEUCHOS_LIBRARY_STATIC)
	set(SACADO_LIBRARIES_STATIC ${SACADO_LIBRARY_STATIC} ${TEUCHOS_LIBRARY_STATIC})
else()
    set(SACADO_LIBRARIES_STATIC ${SACADO_LIBRARY_STATIC})
endif()

if(TEUCHOS_LIBRARY)
	set(SACADO_LIBRARIES ${SACADO_LIBRARY} ${TEUCHOS_LIBRARY})	
else()
	set(SACADO_LIBRARIES ${SACADO_LIBRARY})
endif()
			 		 
if(TEUCHOS_INCLUDE_DIR)
	set(SACADO_INCLUDE_DIRS ${SACADO_INCLUDE_DIR} ${TEUCHOS_INCLUDE_DIR})
else()
	set(SACADO_INCLUDE_DIRS ${SACADO_INCLUDE_DIR})
endif()
	
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sacado  DEFAULT_MSG
	                              SACADO_LIBRARY SACADO_INCLUDE_DIR)
mark_as_advanced(SACADO_INCLUDE_DIR SACADO_LIBRARY )

	
