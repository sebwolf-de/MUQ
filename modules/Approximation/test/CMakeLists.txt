CreateTestGroup(
  APPROXIMATION_GP_Kernels

  GaussianProcesses/KernelTests.cpp
  #GaussianProcesses/KarhunenLoeveTest.cpp
)

CreateTestGroup(
  APPROXIMATION_GP

  GaussianProcesses/GaussianProcessTests.cpp
  GaussianProcesses/StateSpaceTests.cpp
)

CreateTestGroup(
  APPROXIMATION_POLYNOMIALS

  Polynomials/BasisExpansionTest.cpp
  Polynomials/MonotoneExpansionTest.cpp
  Polynomials/PolynomialTests.cpp
)

#CreateTestGroup(
#	APPROXIMATION_REGRESSION
#
#	RegressionTests/PolynomialTests.cpp
#	RegressionTests/MultiIndexTests.cpp
#	RegressionTests/RegressionTests.cpp
#	RegressionTests/LocalRegressionTests.cpp
#)
