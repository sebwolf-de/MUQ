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
    HDF5/BlockDataset.cpp
)

CreateCompileGroup(
    UTILITIES_LINEARALGEBRA
    "Extensions and generalizations of Eigen for linear algebra."
    muqUtilities
    ""
    "EIGEN3"
    "SUNDIALS"
    RandomGenerator.cpp
    LinearAlgebra/LinearOperator.cpp
    LinearAlgebra/CompanionMatrix.cpp
    LinearAlgebra/BlockDiagonalOperator.cpp
    LinearAlgebra/BlockRowOperator.cpp
    LinearAlgebra/KroneckerProductOperator.cpp
    LinearAlgebra/SumOperator.cpp
    LinearAlgebra/ConcatenateOperator.cpp
    LinearAlgebra/DiagonalOperator.cpp
    LinearAlgebra/ProductOperator.cpp

    LinearAlgebra/AnyAlgebra.cpp
    LinearAlgebra/ScalarAlgebra.cpp
    LinearAlgebra/EigenVectorAlgebra.cpp
    LinearAlgebra/EigenMatrixAlgebra.cpp
    LinearAlgebra/SundialsAlgebra.cpp
)
