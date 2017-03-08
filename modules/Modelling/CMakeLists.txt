add_subdirectory(src)


add_subdirectory(test)
set(all_gtest_sources ${all_gtest_sources} ${Test_Modelling_Sources} PARENT_SCOPE)

