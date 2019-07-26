# make sure that the NANOFLANN library is available
set(CMAKE_REQUIRED_LIBRARIES ${NANOFLANN_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${NANOFLANN_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES(
  "
#include <nanoflann.hpp>

#include <stdio.h>

using namespace nanoflann;

int main(int argc, char** argv)
{
  const size_t N = 1000;

  PointCloud<double> cloud;
  generateRandomPointCloud(cloud, N);
  typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, PointCloud<double>>, PointCloud<double>, 3> my_kd_tree_t;

  my_kd_tree_t   index(3, cloud, KDTreeSingleIndexAdaptorParams(10) );
  index.buildIndex();

  return 0;
}

  "
  NANOFLANN_COMPILES)


	if(NOT NANOFLANN_COMPILES)
		set(NANOFLANN_TEST_FAIL 1)
	else()
		set(NANOFLANN_TEST_FAIL 0)
	endif()
