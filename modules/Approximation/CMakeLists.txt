add_subdirectory(src)

if(BUILD_APPROXIMATION_TESTS)

  #pass the test sources to the higher level
  add_subdirectory(test)
  set(all_gtest_sources ${all_gtest_sources} ${APPROXIMATION_TEST_SOURCES} PARENT_SCOPE)

endif()
