########################################
##### LOOK FOR GTEST              ######
########################################
IF(MUQ_USE_GTEST)
	find_package(GTEST)
	
	IF(GTEST_FOUND)
		include(CheckGTEST)	
	Endif()
	
	IF(GTEST_FOUND AND NOT GTEST_TEST_FAIL)

		set(MUQ_BUILD_TESTS ON)
		
		include_directories(${GTEST_INCLUDE_DIRS})
		LIST(APPEND test_link_libs ${GTEST_LIBRARIES})

	else()

		message(WARNING "Could not find GTEST.  No tests can be compiled!")
	    	set(MUQ_BUILD_TESTS OFF)
	
	

	ENDIF(GTEST_FOUND AND NOT GTEST_TEST_FAIL)

else(MUQ_USE_GTEST)
    
    message(WARNING “Tried to compile tests, but MUQ_USE_GTEST is OFF.  Turning off tests.”)
    set(MUQ_BUILD_TESTS OFF)
		
endif(MUQ_USE_GTEST)        

########################################
##### REMOVE DUPLICATE INCLUDES   ######
########################################

list( REMOVE_DUPLICATES MUQ_LINK_LIBS)
