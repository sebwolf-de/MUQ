CreateTestGroup(
    MODELING_CORE

    WorkPieceTests.cpp
    #WorkPieceDerivativesTests.cpp
    WorkGraphTests.cpp
    WorkGraphPieceTests.cpp
    ModPieceTests.cpp
    ConstantPieceTests.cpp
    IdentityPieceTests.cpp

    FlannTest/FlannCacheTests.cpp
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
    DistributionsTests/UniformDistributionTests.cpp
    DistributionsTests/GaussianDistributionTests.cpp
)
