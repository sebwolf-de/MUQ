
########################################
##### LOOK FOR GLOG               ######
########################################
IF(MUQ_USE_GLOG)
	find_package(GLOG)
	
	IF(GLOG_FOUND)
		include(CheckGLOG)
		if(GLOG_TEST_FAIL)
			set(MUQ_USE_GLOG OFF)
		else()
			ADD_DEFINITIONS(-DMUQ_GLOG)
			include_directories(${GLOG_INCLUDE_DIRS})
			LIST(APPEND MUQ_EXTERNAL_INCLUDES ${GLOG_INCLUDE_DIRS})
			LIST(APPEND MUQ_LINK_LIBS ${GLOG_LIBRARIES})
			LIST(APPEND MUQ_LINK_LIBS_STATIC ${GLOG_LIBRARIES_STATIC})
		endif()
	ELSE()
		set(MUQ_USE_GLOG OFF)
		
		# send a warning that no logging will be available
		message(WARNING "Could not find GLOG.  This means that no logging will be available in MUQ!")
		
	ENDIF()

endif(MUQ_USE_GLOG)

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
		include_directories(${GTEST_INCLUDE_DIRS})
		LIST(APPEND test_link_libs ${GTEST_LIBRARIES})
	else()
		# if we couldn't find gtest, turn off all our tests
		    set(Inference_tests OFF)
		    set(Optimization_tests OFF)
		    set(IncrementalApproximation_tests OFF)
		    set(Geostats_tests OFF)
		    set(UtilitiesAndModelling_tests OFF)
		    set(PolynomialChaos_tests OFF)
		    set(Regression_tests OFF)
		    set(PDE_tests OFF)
		
		# send a warning that no tests will be compiled
		message(WARNING "Could not find GTEST.  No tests can be compiled!")

	ENDIF(GTEST_FOUND AND NOT GTEST_TEST_FAIL)

    else(MUQ_USE_GTEST)
        message(WARNING “Tried to compile tests, but MUQ_USE_GTEST is OFF.  Turning off tests.”)
	    set(Inference_tests OFF)
	    set(Optimization_tests OFF)
	    set(IncrementalApproximation_tests OFF)
	    set(Geostats_tests OFF)
	    set(UtilitiesAndModelling_tests OFF)
	    set(PolynomialChaos_tests OFF)
	    set(Regression_tests OFF)
	    set(PDE_tests OFF)
		
    endif(MUQ_USE_GTEST)
endif()


########################################
##### LOOK FOR SACADO             ######
########################################
SET(MUQ_Sacado 0)

if(MUQ_USE_SACADO)
	
    FIND_PACKAGE(SACADO)

    
    IF (SACADO_FOUND)
      ADD_DEFINITIONS(-DMUQ_Sacado)
      SET(MUQ_Sacado 1)
    
      # include the sacado library for linking
      LIST(APPEND MUQ_LINK_LIBS ${SACADO_LIBRARIES})
      LIST(APPEND MUQ_LINK_LIBS_STATIC ${SACADO_LIBRARIES_STATIC})
      include_directories(${SACADO_INCLUDE_DIRS})
      LIST(APPEND MUQ_EXTERNAL_INCLUDES ${SACADO_INCLUDE_DIRS})
      
    ELSE()
      if(PDE_build)
        message("SACADO NOT found, turning off MUQ PDE")
      endif()
      set(MUQ_USE_SACADO OFF)
      set(PDE_build OFF)
      set(PDE_tests OFF)
    ENDIF()
else(MUQ_USE_SACADO)
  set(PDE_build OFF)
  set(PDE_tests OFF)
endif(MUQ_USE_SACADO)

########################################
##### LOOK FOR LIBMESH            ######
########################################
SET(MUQ_LIBMESH 0)

if( MUQ_USE_LIBMESH )
  FIND_PACKAGE(LIBMESH)
  IF( LIBMESH_FOUND )
    ADD_DEFINITIONS(-DMUQ_LIBMESH)

    SET(MUQ_LIBMESH 1)

    # include the sacado library for linking
    LIST(APPEND MUQ_LINK_LIBS ${LIBMESH_LIBRARIES})
    LIST(APPEND MUQ_LINK_LIBS_STATIC ${LIBMESH_LIBRARIES_STATIC})
    include_directories(${LIBMESH_INCLUDE_DIRS})
    LIST(APPEND MUQ_EXTERNAL_INCLUDES ${LIBMESH_INCLUDE_DIRS})
    
  ELSE()
    set(MUQ_USE_LIBMESH OFF)
    set(PDE_build OFF)
    set(PDE_tests OFF)
  ENDIF()
else(MUQ_USE_LIBMESH)
  set(PDE_build OFF)
  set(PDE_tests OFF)
endif( MUQ_USE_LIBMESH )

########################################
##### LOOK FOR CUDA               ######
########################################
if(MUQ_USE_CUDA)
	FIND_PACKAGE(CUDA)
	if(CUDA_FOUND)
		set(CUDA_PROPAGATE_HOST_FLAGS OFF)
		set(CUDA_NVCC_FLAGS "-arch=sm_20")
	else(CUDA_FOUND)
		set(MUQ_USE_CUDA OFF)
	endif(CUDA_FOUND)
endif(MUQ_USE_CUDA)

########################################
##### LOOK FOR ARMADILLO          ######
########################################
if(MUQ_USE_ARMADILLO)
  FIND_PACKAGE(Armadillo)
  if(ARMADILLO_FOUND)
    add_definitions( -DMUQ_Armadillo)
	set(use_armadillo ON)
	
	include_directories(${ARMADILLO_INCLUDE_DIR})
	LIST(APPEND MUQ_EXTERNAL_INCLUDES ${ARMADILLO_INCLUDE_DIR})
	LIST(APPEND MUQ_LINK_LIBS ${ARMADILLO_LIBRARY})
	
  else(ARMADILLO_FOUND)
	set(use_armadillo OFF)
  endif(ARMADILLO_FOUND)
endif(MUQ_USE_ARMADILLO)

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
    
    #FIND_PACKAGE(MKL)
	
    #IF (MKL_FOUND)
		
		# include the mkl library for linking
		if(MUQ_USE_OPENMP)
		  LIST(APPEND MUQ_LINK_LIBS ${MUQ_MKL_DIR}/lib/intel64/libmkl_intel_lp64${CMAKE_SHARED_LIBRARY_SUFFIX} ${MUQ_MKL_DIR}/lib/intel64/libmkl_core${CMAKE_SHARED_LIBRARY_SUFFIX} ${MUQ_MKL_DIR}/lib/intel64/libmkl_gnu_thread${CMAKE_SHARED_LIBRARY_SUFFIX})
		else()
		  LIST(APPEND MUQ_LINK_LIBS ${MUQ_MKL_DIR}/lib/intel64/libmkl_intel_lp64${CMAKE_SHARED_LIBRARY_SUFFIX} ${MUQ_MKL_DIR}/lib/intel64/libmkl_core${CMAKE_SHARED_LIBRARY_SUFFIX} ${MUQ_MKL_DIR}/lib/intel64/libmkl_sequential${CMAKE_SHARED_LIBRARY_SUFFIX})
		endif()
		
		include_directories(${MUQ_MKL_DIR}/include)
		LIST(APPEND MUQ_EXTERNAL_INCLUDES ${MUQ_MKL_DIR}/include)
		add_definitions(-DEIGEN_USE_MKL_ALL)
    #ELSE()
    #message(WARNING "Unable to find Intel MKL")
#		set(MUQ_USE_MKL OFF)
    #ENDIF()
endif()

########################################
##### LOOK FOR PYTHON             ######
########################################

if(MUQ_USE_PYTHON)

  if(DEFINED MUQ_PYTHON_VERSION)

    FIND_PACKAGE(PythonLibs ${MUQ_PYTHON_VERSION} EXACT REQUIRED)

    if(NOT PYTHON_LIBRARY)
    message(WARNING "Could not find an exact version match for Python${MUQ_PYTHON_VERSION}.  MUQ will not be compiled with Python support.") 
      set(MUQ_USE_PYTHON OFF)
    else()
      include_directories(${PYTHON_INCLUDE_DIR})
      LIST(APPEND MUQ_EXTERNAL_INCLUDES ${PYTHON_INCLUDE_DIR})
      LIST(APPEND MUQ_LINK_LIBS ${PYTHON_LIBRARY})
      LIST(APPEND MUQ_LINK_LIBS_STATIC ${PYTHON_LIBRARY_STATIC})
    endif()
    
  else()
    FIND_PACKAGE(PythonLibs 2.7)
    
    if(NOT PYTHON_LIBRARY)
    message(WARNING "Could not find a Python Library.  MUQ will not be compiled with Python support.") 
      set(MUQ_USE_PYTHON OFF)
    else()
      message("PYTHON_INCLUDE_DIR = ${PYTHON_INCLUDE_DIR}")
      include_directories(${PYTHON_INCLUDE_DIR})
      LIST(APPEND MUQ_EXTERNAL_INCLUDES ${PYTHON_INCLUDE_DIR})
      LIST(APPEND MUQ_LINK_LIBS ${PYTHON_LIBRARY})
      LIST(APPEND MUQ_LINK_LIBS_STATIC ${PYTHON_LIBRARY_STATIC})
    endif()
   
  endif()

endif()

########################################
##### REMOVE DUPLICATE INCLUDES   ######
########################################
list( REMOVE_DUPLICATES MUQ_LINK_LIBS)
