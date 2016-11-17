# make sure that we can compile and link to glog
set(CMAKE_REQUIRED_LIBRARIES ${GLOG_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${GLOG_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <glog/logging.h>
  int main(int argc, char* argv[]) {
       google::InitGoogleLogging(argv[0]);
       LOG(INFO) << 10;
     }
  "
  GLOG_COMPILES)
  
  if(NOT GLOG_COMPILES)
  	set(GLOG_TEST_FAIL 1)
  else()
  	set(GLOG_TEST_FAIL 0)
  endif()