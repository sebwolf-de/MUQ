########################################
##### LOOK FOR GTEST              ######
########################################
IF(UtilitiesAndModelling_tests)
    IF(MUQ_USE_GTEST)
         find_package(GTEST)
  
         IF(GTEST_FOUND)
	     include(CheckGTEST)	
	 Endif()
	
	 IF(GTEST_FOUND AND NOT GTEST_TEST_FAIL)
             include_directories(${GTEST_INCLUDE_DIRS})
	     LIST(APPEND test_link_libs ${GTEST_LIBRARIES})
         else()
	     # if we couldn't find gtest, turn off all our tests
	     set(UtilitiesAndModelling_tests OFF)
	     # send a warning that no tests will be compiled
	     message(WARNING "Could not find GTEST.  No tests can be compiled!")
	 ENDIF(GTEST_FOUND AND NOT GTEST_TEST_FAIL)

    else(MUQ_USE_GTEST)
        message(WARNING “Tried to compile tests, but MUQ_USE_GTEST is OFF.  Turning off tests.”)
	set(UtilitiesAndModelling_tests OFF)
    endif(MUQ_USE_GTEST)
endif()