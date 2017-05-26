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
    MODELING_CORE
    "The basic building blocks of the Modeling library."
    muqModeling
    ""
    "EIGEN3;BOOST"
    ""
    WorkPiece.cpp

    NodeNameFinder.cpp
    AnyAlgebra.cpp

    WorkGraph.cpp
    WorkGraphNode.cpp
    WorkGraphEdge.cpp

    WorkGraphPiece.cpp

    ConstantPiece.cpp
    IdentityPiece.cpp
)

CreateCompileGroup(
    SUNDIALS_MODELS
    "ODE integrators using Sundials"
    muqModeling
    "MODELING_CORE"
    "EIGEN3;BOOST;SUNDIALS"
    ""
    ODEBase.cpp
    ODE.cpp
    RootfindingIVP.cpp
)
