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

    LinearAlgebra/ScalarAlgebraTests.cpp
    LinearAlgebra/EigenVectorAlgebraTests.cpp
    LinearAlgebra/EigenMatrixAlgebraTests.cpp
    LinearAlgebra/SundialsAlgebraTests.cpp

    MultiIndices/MultiIndexTests.cpp
    MultiIndices/MultiIndexLimiterTests.cpp
)

CreateTestGroup(
    UTILITIES_MULTIINDEX

    MultiIndices/MultiIndexTests.cpp
)
