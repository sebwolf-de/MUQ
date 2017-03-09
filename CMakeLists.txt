cmake_minimum_required (VERSION 2.8) 

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/finds)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/checks)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/builds)
  
INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckCXXSourceRuns)
INCLUDE(CheckFunctionExists)

# define the project name
project (MUQ) 
set(MUQ_MAJOR_VERSION 1)
set(MUQ_MINOR_VERSION 0)
set(MUQ_PATCH_VERSION 0)
set(MUQ_VERSION ${MUQ_MAJOR_VERSION}.${MUQ_MINOR_VERSION}.${MUQ_PATCH_VERSION})


message(STATUS "\n\nThis is MUQ, the MIT UQ Library!\n\n")

###############################################
# Compiler configuration
###############################################
INCLUDE(SetupCompiler)
INCLUDE(CompilerCheck)
include_directories(.)

###############################################
# User options
###############################################
INCLUDE(OptionsSetup)

###############################################
# Find source files and calculate required dependencies
###############################################
INCLUDE(GetPossibleDependencyList)
INCLUDE(CompileGroup)

message(STATUS "All possible MUQ dependencies = ${POSSIBLE_MUQ_DEPENDENCIES}")

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
    message("SOURCES for ${target} are ${${target}_SOURCES}")
endforeach()


###############################################
# Manage dependencies
###############################################
INCLUDE(DependencySearch)

include_directories(. ./external/include)

# write a configure file containing information about what libraries are available
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/MUQ/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/MUQ/config.h)

###############################################
# set module dependencies
###############################################
INCLUDE(DependencyCheck)

###############################################
# Create BUILD targets for cmake
###############################################
if(MUQ_USE_GTEST)
    message(STATUS "GTEST looks good, so we'll build tests.")
    set(MUQ_BUILD_TESTS ON)
endif()

# Build all the targets
foreach(libName ${MUQ_TARGETS})

    message("LOOKING AT ${libName}")
    message("${${libName}_SOURCES}")
    
    list(LENGTH ${libName}_SOURCES strLength)
    if(${strLength} GREATER 0)
        message("ADDING TARGET ${libName}")
	message("ADDING SOURCES ${${libName}_SOURCES}")
	
        ADD_LIBRARY(${libName} ${${libName}_SOURCES})
        TARGET_LINK_LIBRARIES(${libName} ${MUQ_LINK_LIBS})

        install(TARGETS ${libName}
                EXPORT MUQDepends
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



# Create the list of libraries we just built
INCLUDE(BuiltLibrarySetup)

###############################################################
# INSTALLATION INFORMATION
###############################################################
get_directory_property(MUQ_COMPILE_DEFINITIONS COMPILE_DEFINITIONS)
message(STATUS "Storing MUQ_COMPILE_DEFINITIONS = ${MUQ_COMPILE_DEFINITIONS}")

set(MUQ_EXPORT_LIBRARIES "")
FOREACH(BUILT_LIB ${MUQ_BUILT_LIBRARIES})
	list(APPEND MUQ_EXPORT_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/${BUILT_LIB})
ENDFOREACH()

export(PACKAGE MUQ)

# Create the MUQConfig.cmake file
configure_file(cmake/MUQConfig.cmake.in "${PROJECT_BINARY_DIR}/MUQConfig.cmake" @ONLY)
configure_file(cmake/MUQConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/MUQConfig.cmake" @ONLY)

# Install the MUQConfig.cmake file
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/MUQConfig.cmake" 
        DESTINATION "${CMAKE_INSTALL_PREFIX}/CMake/MUQ"
	    COMPONENT dev)
  
# Install the export set for use with the install-tree
install(EXPORT MUQDepends
        DESTINATION "${CMAKE_INSTALL_PREFIX}/CMake/MUQ"
	    COMPONENT dev)
  
# install all the MUQ headers -- may be better off installing for each package individually
install(DIRECTORY MUQ
        DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
	    FILES_MATCHING PATTERN "*.h")
  
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${INSTALL_LIB_DIR}) 
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${INSTALL_LIB_DIR}) 



###############################################
# Write summary to terminal
###############################################

# write a summary of the build configuration
include(MakeSummaryFile)

FILE(READ ${CMAKE_BINARY_DIR}/summary.log MUQ_LOG_SUMMARY)
MESSAGE("${MUQ_LOG_SUMMARY}")


###############################################
# Add a target to build documentation
###############################################
# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
	message(STATUS "To build documentation, run 'make doc'")
	
	option(MUQDOC_INCLUDE_ANALYTICS OFF)
	if(MUQDOC_INCLUDE_ANALYTICS)
		message(STATUS "Will include google analytics header information in doxygen output.")
		set(ANALYTICS_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/MUQdocumentation/doxFiles/analytics_header.html)
	endif()
	  
  message(STATUS "Extracting parameter values from cpp files.  Results are stored in ${CMAKE_BINARY_DIR}/ExtractedParameters.dox")
  message(STATUS "python ${CMAKE_CURRENT_SOURCE_DIR}/MUQdocumentation/doxFiles/FindPtreeParameters.py ${CMAKE_CURRENT_SOURCE_DIR}")
  execute_process(COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/MUQdocumentation/doxFiles/FindPtreeParameters.py ${CMAKE_CURRENT_SOURCE_DIR}
                  OUTPUT_FILE "${CMAKE_BINARY_DIR}/ParameterDescriptions.dox")
    
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/MUQdocumentation/muq.doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/muq.doxyfile @ONLY)
    #configure_file(${CMAKE_CURRENT_SOURCE_DIR}/MUQdocumentation/ParameterDescriptions.dox.in ${CMAKE_CURRENT_BINARY_DIR}/ParameterDescriptions.dox @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/muq.doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
    
    
endif(DOXYGEN_FOUND)

