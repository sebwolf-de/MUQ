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

    # Make sure all upstream dependency groups are enabled  
    message(STATUS "Configuring compile group ${group}")
    foreach(depend ${${group}_REQUIRES_GROUPS})
        if(MUQ_ENABLEGROUP_${group} AND NOT MUQ_ENABLEGROUP_${depend})
            message(STATUS "    The ${group} group depends on the ${depend} group, but the ${depend} group was not enabled.")
            message(STATUS "    Turning the ${depend} group on.")
            set(MUQ_ENABLEGROUP_${depend} ON)
        endif()
    endforeach()

    # Add to the list of required external libraries
    foreach(depend ${${group}_REQUIRES})
        list(APPEND MUQ_REQUIRES ${depend})
    endforeach()

    # Add to the list of desired (i.e., optional) external libraries
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
        if(MUQ_ENABLEGROUP_${group})
	    if(${${group}_LIBRARY} MATCHES ${target})

                # Check to see if a group has any source (e.g., *.cpp) files.  Flag it as something that will be built if it does.
	        list(LENGTH ${group}_SOURCES sources_length)
		if(sources_length GREATER 0)
	            set(${group}_IS_COMPILED ON CACHE INTERNAL "Whether or not the ${group} is used in any library.")
		endif()
		
	        list(APPEND ${target}_SOURCES ${${group}_SOURCES})
	    endif()
	endif()
    endforeach()

    if(${target}_SOURCES)
        list(REMOVE_DUPLICATES ${target}_SOURCES)
    endif()
    
endforeach()
