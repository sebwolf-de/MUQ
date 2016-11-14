# begin accumulating the tests
set(all_gtest_sources RunTests.cpp)

set(test_link_libs )

if(UtilitiesAndModelling_build)
    add_subdirectory(Modelling)

    # always link in modelling and general utilities
    list(APPEND test_link_libs muqModelling)
endif(UtilitiesAndModelling_build)

#only build the tests if some of them should be built
IF(UtilitiesAndModelling_tests)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
  
    CHECK_CXX_COMPILER_FLAG("-std=c++11" HAS_PTHREAD)
    if(HAS_PTHREAD)
	set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pthread")
    endif()
  
    ADD_EXECUTABLE(RunAllTests ${all_gtest_sources})
    message("MUQ_LINK_LIBS=${MUQ_LINK_LIBS}")
    TARGET_LINK_LIBRARIES(RunAllTests ${test_link_libs} ${GTEST_LIBRARY} ${external_test_libs})
ENDIF()
