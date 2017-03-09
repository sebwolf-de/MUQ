function(CreateCompileGroup GROUP_NAME LIBRARY_NAME DEPENDENCIES)

  set(MUQ_GROUPS "${MUQ_GROUPS};${GROUP_NAME}" CACHE INTERNAL "A list of all of the muq groups.")
  set(${GROUP_NAME}_REQUIRES ${DEPENDENCIES} CACHE INTERNAL "External packages required by the ${GROUP_NAME} compile group.")


  # Compute the path to the source file relative to the root directory
  string(REPLACE "${CMAKE_SOURCE_DIR}/"  ""  RELATIVE_DIR ${CMAKE_CURRENT_LIST_DIR})

  set(SOURCES )
  foreach(source ${ARGN})
    list(APPEND SOURCES "${RELATIVE_DIR}/${source}")
  endforeach()

  set(${GROUP_NAME}_SOURCES ${SOURCES} CACHE INTERNAL "Source files included in the ${GROUP_NAME} compile group.")

endfunction(CreateCompileGroup)