add_subdirectory(src)

if(MUQ_BUILD_TESTS)
    add_subdirectory(test)
    set(all_gtest_sources ${all_gtest_sources} ${Test_Modelling_Sources} PARENT_SCOPE)
endif()
