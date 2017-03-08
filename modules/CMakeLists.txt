# begin accumulating the tests
set(all_gtest_sources RunTests.cpp)

set(test_link_libs )

add_subdirectory(Modelling)
add_subdirectory(Approximation)


list(APPEND test_link_libs muqModelingCore)


#only build the tests if some of them should be built
IF(MUQ_BUILD_TESTS)

    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
  
    CHECK_CXX_COMPILER_FLAG("-std=c++11" HAS_PTHREAD)
    if(HAS_PTHREAD)
	set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pthread")
    endif()
  
    ADD_EXECUTABLE(RunAllTests ${all_gtest_sources})

    message("MUQ_LINK_LIBS=${MUQ_LINK_LIBS}")
    TARGET_LINK_LIBRARIES(RunAllTests ${test_link_libs} ${GTEST_LIBRARY} ${external_test_libs})

ENDIF()
