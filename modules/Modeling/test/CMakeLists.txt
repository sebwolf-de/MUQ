CreateTestGroup(
    MODELING_CORE
    
    WorkPieceTests.cpp
    WorkPieceDerivativesTests.cpp
    WorkGraphTests.cpp
    WorkGraphPieceTests.cpp

    ConstantPieceTests.cpp
    IdentityPieceTests.cpp
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

