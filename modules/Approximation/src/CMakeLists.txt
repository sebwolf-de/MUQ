# USAGE:
# CreateCompileGroup(
#     <group name>
#     <brief description>
#     <library_name>
#     <other compile group dependencies>
#     <required dependencies>
#     <optional dependencies>
#     <source file 1>
#     <source file 2>
#     ...
#     <source file N>
# )

CreateCompileGroup(
	APPROXIMATION_REGRESSION
	"Polynomial regression (global and local)"
	muqApproximation
	"MODELING_FLANN;MODELING_CORE;APPROXIMATION_POLYNOMIALS"
	""
	""
	Regression/MultiIndex.cpp
	
	Regression/Regression.cpp
	Regression/LocalRegression.cpp
)


CreateCompileGroup(
	APPROXIMATION_POLYNOMIALS
	"Polynomials and multiindices"
	muqApproximation
	"MODELING_CORE"
	""
	""
	Polynomials/Polynomial.cpp
	Polynomials/Monomial.cpp
	Polynomials/PhysicistHermite.cpp
        Polynomials/ProbabilistHermite.cpp
	Polynomials/Legendre.cpp
        Polynomials/Laguerre.cpp
)

CreateCompileGroup(
    APPROXIMATION_GP_Kernels
    "Gaussian process covariance kernels."
    muqApproximation
    "MODELING_CORE"
    "EIGEN3;BOOST"
    ""
)

CreateCompileGroup(
    APPROXIMATION_GP
    "Tools for constructing Gaussian processes."
    muqApproximation
    "MODELING_CORE;APPROXIMATION_GP_Kernels"
    "EIGEN3;BOOST"
    "NLOPT"
)
