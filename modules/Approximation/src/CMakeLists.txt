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
    APPROXIMATION_GP
    "Tools for constructing Gaussian processes."
    muqApproximation
    "MODELING_CORE"
    "EIGEN3"
    "NLOPT"
    GaussianProcesses/GaussianProcess.cpp
)
