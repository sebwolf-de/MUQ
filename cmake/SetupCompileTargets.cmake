# PURPOSE:
# This file sets up the MUQ build targets (e.g., libMuqModeling).  Information is
# used from the compile groups that were processed in the ProcessCompileGroups.cmake
# file.
#

set(MUQ_LIBRARIES )

# Build all the targets
foreach(libName ${MUQ_TARGETS})

    list(LENGTH ${libName}_SOURCES strLength)
    if(${strLength} GREATER 0)
    
        ADD_LIBRARY(${libName} ${${libName}_SOURCES})
        TARGET_LINK_LIBRARIES(${libName} ${${CMAKE_PROJECT_NAME}_LINK_LIBS})

        # Add dependencies for any required dependencies that MUQ is going to build internally
        foreach(depend ${MUQ_REQUIRES})
            message(STATUS "Checking for dependency of ${libName} on internal build of ${depend}")
            if(USE_INTERNAL_${depend})
                message(STATUS "Adding dependency of ${libName} on ${depend}")
                add_dependencies(${libName} ${depend})
            endif()
        endforeach()
    
        list(APPEND MUQ_LIBRARIES ${libName})
        install(TARGETS ${libName}
                EXPORT ${CMAKE_PROJECT_NAME}Depends
                LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
                ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
    endif()
    
endforeach()

# If a group depends on an external library that is going to be built by MUQ, then make sure we account for that dependency
foreach(group ${MUQ_GROUPS})

    list(LENGTH ${group}_SOURCES strLength)

    foreach(depend ${POSSIBLE_MUQ_DEPENDENCIES})
        list(FIND ${group}_REQUIRES ${depend} needsExternal)

        if(USE_INTERNAL_${depend})
            if(needsExternal AND ${USE_INTERNAL_${depend}} AND (strLength GREATER 0))
                add_dependencies(${${group}_LIBRARY} ${depend})
            endif()
	endif()
    endforeach()
endforeach()
