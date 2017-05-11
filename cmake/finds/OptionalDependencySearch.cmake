
########################################
##### LOOK FOR GTEST              ######
########################################
IF(Inference_tests OR Optimization_tests OR IncrementalApproximation_tests OR Geostats_tests 
  OR UtilitiesAndModelling_tests OR polychaos_tests OR Regression_tests)
  
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
endif()

########################################
##### LOOK FOR NLOPT              ######
########################################
if(MUQ_USE_NLOPT)
        FIND_PACKAGE(NLOPT)
	
	IF (NLOPT_FOUND)
		add_definitions(-DMUQ_USE_NLOPT)
		# include the sacado library for linking
		LIST(APPEND MUQ_LINK_LIBS ${NLOPT_LIBRARIES})
		LIST(APPEND MUQ_LINK_LIBS_STATIC ${NLOPT_LIBRARIES_STATIC})
		
		include_directories(${NLOPT_INCLUDE_DIRS})
		LIST(APPEND MUQ_EXTERNAL_INCLUDES ${NLOPT_INCLUDE_DIRS})
		
    ELSE()
		set(MUQ_USE_NLOPT OFF)
    ENDIF()

endif(MUQ_USE_NLOPT)


########################################
##### LOOK FOR MKL                ######
########################################
if(MUQ_USE_MKL)
    # include the mkl library for linking
    if(MUQ_USE_OPENMP)
	  LIST(APPEND MUQ_LINK_LIBS ${MUQ_MKL_DIR}/lib/intel64/libmkl_intel_lp64${CMAKE_SHARED_LIBRARY_SUFFIX} ${MUQ_MKL_DIR}/lib/intel64/libmkl_core${CMAKE_SHARED_LIBRARY_SUFFIX} ${MUQ_MKL_DIR}/lib/intel64/libmkl_gnu_thread${CMAKE_SHARED_LIBRARY_SUFFIX})
    else()
	  LIST(APPEND MUQ_LINK_LIBS ${MUQ_MKL_DIR}/lib/intel64/libmkl_intel_lp64${CMAKE_SHARED_LIBRARY_SUFFIX} ${MUQ_MKL_DIR}/lib/intel64/libmkl_core${CMAKE_SHARED_LIBRARY_SUFFIX} ${MUQ_MKL_DIR}/lib/intel64/libmkl_sequential${CMAKE_SHARED_LIBRARY_SUFFIX})
    endif()
		
    include_directories(${MUQ_MKL_DIR}/include)
    LIST(APPEND MUQ_EXTERNAL_INCLUDES ${MUQ_MKL_DIR}/include)
    add_definitions(-DEIGEN_USE_MKL_ALL)

endif()

########################################
##### REMOVE DUPLICATE INCLUDES   ######
########################################
list( REMOVE_DUPLICATES MUQ_LINK_LIBS)
