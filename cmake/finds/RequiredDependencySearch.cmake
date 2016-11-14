# define a macro to look for a package and install a local copy if we can't find it
macro(GetDependency name)
    find_package(${name})
    if(${name}_FOUND)
        # check to make sure the library can be linked to
	include(Check${name})

   	if(NOT ${name}_TEST_FAIL)
	    set(USE_INTERNAL_${name} 0)
	 else()
	    set(USE_INTERNAL_${name} 1)	
	 endif()

    else()
        set(USE_INTERNAL_${name} 1) 
    endif()

    if(USE_INTERNAL_${name})
        include(Build${name})
    endif()
													
    # store include directory information
    include_directories(${${name}_INCLUDE_DIRS})
    LIST(APPEND MUQ_EXTERNAL_INCLUDES ${${name}_INCLUDE_DIRS})
    # store library information
    LIST(APPEND MUQ_LINK_LIBS ${${name}_LIBRARIES})
    LIST(APPEND MUQ_LINK_LIBS_STATIC ${${name}_LIBRARIES_STATIC})
endmacro(GetDependency)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/include)

########################################
##### LOOK FOR AND/OR BUILD Eigen ######
########################################
GetDependency(EIGEN3)
