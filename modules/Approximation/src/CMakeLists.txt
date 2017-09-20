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
	
	Regression/Regression.cpp
)
