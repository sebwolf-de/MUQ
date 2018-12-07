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
    OPTIMALEXPERIMENTALDESIGN_CORE_PYTHON
    "Python wrappers for the basic building blocks of the OptimalExperimentalDesign library."
    pymuqOptimalExperimentalDesign
    "OPTIMALEXPERIMENTALDESIGN_CORE"
    "PYTHON"
    ""
    Module.cpp
    EvidenceWrapper.cpp
    OEDResidualWrapper.cpp
    UtilityWrapper.cpp
)
