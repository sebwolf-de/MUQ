CreateTestGroup(
    MODELING_CORE

    WorkPieceTests.cpp
    #WorkPieceDerivativesTests.cpp
    WorkGraphTests.cpp
    WorkGraphPieceTests.cpp
    ModGraphPieceTests.cpp
    ModPieceTests.cpp
    ConstantPieceTests.cpp
    IdentityPieceTests.cpp

    FlannTest/FlannCacheTests.cpp

    CwiseOperatorTests/UnaryOperatorTest.cpp
)

CreateTestGroup(
    MODELING_SUNDIALS_MODELS

    ODETests.cpp
    RootfindingIVPTests.cpp
)

CreateTestGroup(
    MODELING_FLANN

    FlannTest/FlannCacheTests.cpp
)

CreateTestGroup(
    MODELING_DISTRIBUTIONS

    DistributionsTests/DistributionTests.cpp
    DistributionsTests/InverseGammaTests.cpp
    DistributionsTests/UniformDistributionTests.cpp
    DistributionsTests/GaussianDistributionTests.cpp
    DistributionsTests/DensityProductTests.cpp
)
