cmake_minimum_required (VERSION 2.8) 

find_package(MUQ_PRIVATE REQUIRED)

include_directories(${MUQ_PRIVATE_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS ${MUQ_PRIVATE_CXX_FLAGS})
set(CMAKE_CXX_COMPILER ${MUQ_PRIVATE_CXX_COMPILER})

add_executable(example-0 example-0.cpp)

target_link_libraries(example-0 ${MUQ_PRIVATE_LIBRARIES} ${MUQ_PRIVATE_LINK_LIBRARIES})
