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
    OPTIMALEXPERIMENTALDESIGN_CORE
    "Framework for optimal experimental design."
    muqOptimalExperimentalDesign
    "MODELING_DISTRIBUTIONS;SAMPLING_ALGORITHM"
    "EIGEN3"
    ""
    Evidence.cpp
    LogDifference.cpp
    Utility.cpp
    OEDResidual.cpp
)
