#only build the tests if some of them should be built
IF(MUQ_USE_GTEST)

    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
  
    CHECK_CXX_COMPILER_FLAG("-std=c++11" HAS_PTHREAD)
    if(HAS_PTHREAD)
	set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pthread")
    endif()

    # Add all of the relevant GTEST sources
    set(all_gtest_sources modules/RunTests.cpp)
    set(all_compiled_libraries )
    
    foreach(group ${MUQ_TEST_GROUPS})
        message("${group}_TEST_SOURCES = ${${group}_TEST_SOURCES}")
	if(${group}_IS_COMPILED)
	    list(APPEND all_compiled_libraries ${${group}_LIBRARY})
        endif()
        if(${MUQ_ENABLEGROUP_${group}})
            list(APPEND all_gtest_sources ${${group}_TEST_SOURCES})
        endif()
    endforeach()

    list(REMOVE_DUPLICATES all_compiled_libraries)
    list(REMOVE_DUPLICATES all_gtest_sources)

    message("ALL TEST SOURCES = ${all_gtest_sources}")
    ADD_EXECUTABLE(RunAllTests ${all_gtest_sources})

    # Make sure the test executable depends on all of the targets
    foreach(target ${all_compiled_libraries})
        add_dependencies(RunAllTests ${target})
    endforeach()

    message("MUQ_LINK_LIBS = ${MUQ_LINK_LIBS}")
    TARGET_LINK_LIBRARIES(RunAllTests ${MUQ_LIBRARIES} ${MUQ_LINK_LIBS} ${GTEST_LIBRARY})

ENDIF()
