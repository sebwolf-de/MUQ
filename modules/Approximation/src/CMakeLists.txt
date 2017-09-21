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
	REGRESSION
	"Polynomial regression (global and local)"
	muqApproximation
	"MODELING_CORE"
	""
	""
	Regression/Polynomial.cpp
	Regression/Monomial.cpp
	Regression/Hermite.cpp
	Regression/Legendre.cpp
	
	Regression/Regression.cpp
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
