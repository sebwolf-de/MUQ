CreateTestGroup(
    UTILITIES_HDF5

    HDF5/DatasetTests.cpp
    HDF5/PathTests.cpp
    HDF5/HDF5FileTest.cpp
)

CreateTestGroup(
    UTILITIES_LINEARALGEBRA
    
    LinearAlgebra/LinearOperatorTests.cpp
    LinearAlgebra/LyapunovSolverTests.cpp
    RandomGeneratorTest.cpp
    
    LinearAlgebra/AnyAlgebraTests.cpp
)
