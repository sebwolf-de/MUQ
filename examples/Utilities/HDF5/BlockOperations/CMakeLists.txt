cmake_minimum_required(VERSION 2.8)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) # Only look for MUQ when building example as standalone project
  find_package(MUQ REQUIRED)
endif()

set(CMAKE_CXX_COMPILER ${MUQ_CXX_COMPILER})
set(CMAKE_CXX_FLAGS ${MUQ_CXX_FLAGS})

project(SimpleReadWrite)

include_directories(${MUQ_INCLUDE_DIRS})

add_executable(BlockOperations BlockOperations.cpp)
target_link_libraries(BlockOperations ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})
