# PURPOSE:
# The code in this file loops through all of the compile groups, extracts all the required
# and optional dependencies, and collects all of the source files needed for  the compile
# targets (e.g., libmuqModeling, etc...)
#



# Initially, we have no targets to build
set(MUQ_TARGETS "" CACHE INTERNAL "List of MUQ libraries to build.")

# Go compile everything
add_subdirectory(modules)

## Figure out what dependencies we actually need
set(MUQ_REQUIRES )
set(MUQ_DESIRES )
foreach(group ${MUQ_GROUPS})

    foreach(depend ${${group}_REQUIRES})
        list(APPEND MUQ_REQUIRES ${depend})
    endforeach()

    foreach(depend ${${group}_DESIRES})
        list(APPEND MUQ_DESIRES ${depend})
    endforeach()
    
endforeach()

# Remove duplicate requirements
list(REMOVE_DUPLICATES MUQ_REQUIRES)
if(MUQ_DESIRES)
    list(REMOVE_DUPLICATES MUQ_DESIRES)
endif()

# Create a list of all MUQ libraries to build
set(MUQ_TARGETS )
foreach(group ${MUQ_GROUPS})
    list(APPEND MUQ_TARGETS ${${group}_LIBRARY})
endforeach()
list(REMOVE_DUPLICATES MUQ_TARGETS)

# Set up the source for each target library
foreach(target ${MUQ_TARGETS})
    set(${target}_SOURCES )

    foreach(group ${MUQ_GROUPS})
        if(MUQ_GROUP_${group})
	    if(${${group}_LIBRARY} MATCHES ${target})
	        list(APPEND ${target}_SOURCES ${${group}_SOURCES})
	    endif()
	endif()
    endforeach()

    if(${target}_SOURCES)
        list(REMOVE_DUPLICATES ${target}_SOURCES)
    endif()
    
endforeach()
