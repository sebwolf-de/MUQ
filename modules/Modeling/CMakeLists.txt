add_subdirectory(src)

if(UtilitiesAndModelling_tests)
  #pass the test sources to the higher level
  add_subdirectory(test)
  set(all_gtest_sources ${all_gtest_sources} ${Test_Modelling_Sources} PARENT_SCOPE)
endif(UtilitiesAndModelling_tests)
