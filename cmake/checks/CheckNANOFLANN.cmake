# make sure that the NANOFLANN library is available
set(CMAKE_REQUIRED_LIBRARIES ${NANOFLANN_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${NANOFLANN_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES(
  "
#include <nanoflann.hpp>
#include <Eigen/Dense>

#include <stdio.h>

using namespace Eigen;
using namespace nanoflann;

int main(int argc, char **argv) {

  const int nSamples = 1000;
  const int dim = 3;
  Eigen::MatrixXd mat = Eigen::MatrixXd::Random(nSamples, dim);

  typedef KDTreeEigenMatrixAdaptor<Eigen::MatrixXd> my_kd_tree_t;

  my_kd_tree_t mat_index(dim, std::cref(mat), 10 /* max leaf */);
  mat_index.index->buildIndex();

  return 0;
}
  "
  NANOFLANN_COMPILES)


	if(NOT NANOFLANN_COMPILES)
		set(NANOFLANN_TEST_FAIL 1)
	else()
		set(NANOFLANN_TEST_FAIL 0)
	endif()
