
find_package(PkgConfig)

set(BOOST_MIN_VERSION "1.56.0")

if(NOT DEFINED MUQ_BOOST_DIR)
	
	unset(Boost_LIBRARIES)
	unset(Boost_INCLUDE_DIR)
	unset(Boost_LIBRARY_DIRS)

	set(Boost_USE_STATIC_LIBS ON)

	if(MUQ_USE_OPENMPI AND MUQ_USE_PYTHON)
		find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS system filesystem regex serialization graph date_time mpi python)
	elseif(MUQ_USE_OPENMPI AND NOT MUQ_USE_PYTHON)
		find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS system filesystem regex serialization graph date_time mpi)
	elseif(NOT MUQ_USE_OPENMPI AND MUQ_USE_PYTHON)
	  find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS system filesystem regex serialization graph date_time python)
	else()
		find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS system filesystem regex serialization graph date_time) 
	endif()

	IF(Boost_FOUND)
		set(BOOST_LIBRARIES_STATIC ${Boost_LIBRARIES})
	endif()


	unset(Boost_LIBRARIES)
	unset(Boost_INCLUDE_DIR)
	unset(Boost_LIBRARY_DIRS)
	unset(Boost_USE_STATIC_LIBS)

	if(MUQ_USE_OPENMPI AND MUQ_USE_PYTHON)
		find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS system filesystem regex serialization graph date_time mpi python)
	elseif(MUQ_USE_OPENMPI AND NOT MUQ_USE_PYTHON)
		find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS system filesystem regex serialization graph date_time mpi)
	elseif(NOT MUQ_USE_OPENMPI AND MUQ_USE_PYTHON)
		find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS system filesystem regex serialization graph date_time python)
	else()
		find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS system filesystem regex serialization graph date_time) 
	endif()

	IF(Boost_FOUND)
		set(BOOST_LIBRARY ${Boost_LIBRARIES})
		set(BOOST_INCLUDE_DIR ${Boost_INCLUDE_DIR})
	endif()

else()

	find_path(BOOST_INCLUDE_DIR boost/property_tree/ptree.hpp
	          HINTS ${MUQ_BOOST_DIR}/include ${MUQ_BOOST_DIR}
	          PATH_SUFFIXES boost NO_DEFAULT_PATH)

	find_library(BOOST_SYSTEM_LIBRARY_STATIC NAMES ${library_prefix}boost_system.${static_library_suffix}
	             HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
	find_library(BOOST_FILESYSTEM_LIBRARY_STATIC NAMES ${library_prefix}boost_filesystem.${static_library_suffix}
	             HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
	find_library(BOOST_TIMER_LIBRARY_STATIC NAMES ${library_prefix}boost_timer.${static_library_suffix}
			 	 HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
	find_library(BOOST_CHRONO_LIBRARY_STATIC NAMES ${library_prefix}boost_chrono.${static_library_suffix}
			 	 HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
 	find_library(BOOST_REGEX_LIBRARY_STATIC NAMES ${library_prefix}boost_regex.${static_library_suffix}
 	             HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
 	find_library(BOOST_SERIAL_LIBRARY_STATIC NAMES ${library_prefix}boost_serialization.${static_library_suffix}
 	             HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
	find_library(BOOST_WSERIAL_LIBRARY_STATIC NAMES ${library_prefix}boost_wserialization.${static_library_suffix}
			  	 HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
 	find_library(BOOST_GRAPH_LIBRARY_STATIC NAMES ${library_prefix}boost_graph.${static_library_suffix}
 	             HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
  	find_library(BOOST_MATH_LIBRARY_STATIC NAMES ${library_prefix}boost_math_tr1.${static_library_suffix}
  	             HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
	find_library(BOOST_MATHF_LIBRARY_STATIC NAMES ${library_prefix}boost_math_tr1f.${static_library_suffix}
			   	 HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
	find_library(BOOST_MATHL_LIBRARY_STATIC NAMES ${library_prefix}boost_math_tr1l.${static_library_suffix}
				 HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
				 
 	find_library(BOOST_SYSTEM_LIBRARY NAMES boost_system
 		 	     HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)	
 	find_library(BOOST_FILESYSTEM_LIBRARY NAMES boost_filesystem
 		 	     HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)	
	find_library(BOOST_TIMER_LIBRARY NAMES boost_timer
			  	 HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
	find_library(BOOST_CHRONO_LIBRARY NAMES boost_chrono
			  	 HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)		
	find_library(BOOST_REGEX_LIBRARY NAMES boost_regex
		 	     HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)	
 	find_library(BOOST_SERIAL_LIBRARY NAMES boost_serialization
 		 	     HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
	find_library(BOOST_WSERIAL_LIBRARY NAMES boost_wserialization
			  	 HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
 	find_library(BOOST_GRAPH_LIBRARY NAMES boost_graph
 		 	     HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)	 
 	find_library(BOOST_MATH_LIBRARY NAMES boost_math_tr1
 		 	     HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
	find_library(BOOST_MATHF_LIBRARY NAMES boost_math_tr1f
			  	 HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
	find_library(BOOST_MATHL_LIBRARY NAMES boost_math_tr1l
				 HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
			 	
	if(MUQ_USE_OPENMPI)
		find_library(BOOST_MPI_LIBRARY NAMES boost_mpi
	        	     HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)	
	   	find_library(BOOST_MPI_LIBRARY_STATIC NAMES ${library_prefix}boost_mpi.${static_library_suffix}
	   	             HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)	 	 
	endif()

	if(MUQ_USE_PYTHON)
		find_library(BOOST_PYTHON_LIBRARY NAMES boost_python
	        	     HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)	
	   	find_library(BOOST_PYTHON_LIBRARY_STATIC NAMES ${library_prefix}boost_python.${static_library_suffix}
	   	             HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)	 	 
	endif()

	set(BOOST_LIBRARY ${BOOST_SYSTEM_LIBRARY} ${BOOST_FILESYSTEM_LIBRARY} ${BOOST_TIMER_LIBRARY} ${BOOST_CHRONO_LIBRARY} ${BOOST_REGEX_LIBRARY} ${BOOST_SERIAL_LIBRARY} ${BOOST_WSERIAL_LIBRARY} ${BOOST_GRAPH_LIBRARY} ${BOOST_MATH_LIBRARY} ${BOOST_MATHF_LIBRARY} ${BOOST_MATHL_LIBRARY} ${BOOST_MPI_LIBRARY} ${BOOST_PYTHON_LIBRARY})
	set(BOOST_LIBRARY_STATIC ${BOOST_SYSTEM_LIBRARY_STATIC} ${BOOST_FILESYSTEM_LIBRARY_STATIC} ${BOOST_TIMER_LIBRARY_STATIC} ${BOOST_CHRONO_LIBRARY_STATIC} ${BOOST_REGEX_LIBRARY_STATIC} ${BOOST_SERIAL_LIBRARY_STATIC} ${BOOST_WSERIAL_LIBRARY_STATIC} ${BOOST_GRAPH_LIBRARY_STATIC} ${BOOST_MATH_LIBRARY_STATIC} ${BOOST_MATHF_LIBRARY_STATIC} ${BOOST_MATHL_LIBRARY_STATIC} ${BOOST_MPI_LIBRARY_STATIC} ${BOOST_PYTHON_LIBRARY})
endif()


set(BOOST_INCLUDE_DIRS ${BOOST_INCLUDE_DIR} )
set(BOOST_LIBRARIES ${BOOST_LIBRARY})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(BOOST  DEFAULT_MSG
                                  BOOST_LIBRARY BOOST_INCLUDE_DIR)

mark_as_advanced(BOOST_INCLUDE_DIR BOOST_LIBRARY )
