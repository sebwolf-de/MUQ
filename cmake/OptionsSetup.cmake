# do you want to build the libraries as shared libraries?
set(BUILD_SHARED_LIBS ON)
option(build_static_libs "Generate additional static libraries" OFF) # generates a couple additional static libs


if (NOT CMAKE_BUILD_TYPE)
 message(STATUS "No build type selected, default to Release")
 set(CMAKE_BUILD_TYPE "Release")
endif()

# set up internal options, such as what modules to build
include(InternalOptions)

# set up external options, such as what optional libraries to use 
include(ExternalOptions)
