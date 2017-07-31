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
    UTILITIES_HDF5
    "MUQ's wrapper fo interacting with HDF5 binary files."
    muqUtilities
    ""
    "EIGEN3;HDF5"
    ""
    HDF5/Attributes.cpp
    HDF5/H5Object.cpp
    HDF5/HDF5File.cpp
    HDF5/HDF5Types.cpp
    HDF5/PathTools.cpp
)
