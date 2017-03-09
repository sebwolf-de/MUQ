function(CreateCompileGroup GROUP_NAME LIBRARY_NAME REQUIRED_DEPENDENCIES OPTIONAL_DEPENDENCIES)

  option(MUQ_GROUP_${GROUP_NAME} "Should the group ${GROUP_NAME} be compiled?" ON)
  
  set(MUQ_GROUPS "${MUQ_GROUPS};${GROUP_NAME}" CACHE INTERNAL "A list of all of the muq groups.")
  
  set(${GROUP_NAME}_REQUIRES ${REQUIRED_DEPENDENCIES} CACHE INTERNAL "External packages required by the ${GROUP_NAME} compile group.")
  set(${GROUP_NAME}_DESIRES ${OPTIONAL_DEPENDENCIES} CACHE INTERNAL "External packages that the ${GROUP_NAME} compile group can exploit if they're available.")
  
  set(${GROUP_NAME}_LIBRARY ${LIBRARY_NAME} CACHE INTERNAL "The library this group will contribute to.")

  # Compute the path to the source file relative to the root directory
  string(REPLACE "${CMAKE_SOURCE_DIR}/"  ""  RELATIVE_DIR ${CMAKE_CURRENT_LIST_DIR})

  set(SOURCES )
  foreach(source ${ARGN})
    list(APPEND SOURCES "${RELATIVE_DIR}/${source}")
  endforeach()

  set(${GROUP_NAME}_SOURCES ${SOURCES} CACHE INTERNAL "Source files included in the ${GROUP_NAME} compile group.")

endfunction(CreateCompileGroup)