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
    "UTILITIES_LINEARALGEBRA"
    "EIGEN3;BOOST;SUNDIALS"
    ""
    WorkPiece.cpp

    NodeNameFinder.cpp

    WorkGraph.cpp
    WorkGraphNode.cpp
    WorkGraphEdge.cpp

    PyModPiece.cpp

    ModPiece.cpp
    ModGraphPiece.cpp
    WorkGraphPiece.cpp
    LinearSDE.cpp

    ConstantPiece.cpp
    ConstantVector.cpp
    IdentityPiece.cpp

    IdentityOperator.cpp
    ReplicateOperator.cpp

    LinearSDE.cpp

    Flann/FlannCache.cpp
)

CreateCompileGroup(
    MODELING_DISTRIBUTIONS
    "Implement probability distributions, which allow us to evaluate log-densities and/or sample the distirubtion."
    muqModeling
    ""
    "MODELING_CORE"
    ""

    Distributions/Distribution.cpp
    Distributions/UniformBox.cpp
    Distributions/Gaussian.cpp
    Distributions/InverseGamma.cpp
    Distributions/DensityProduct.cpp
    Distributions/Density.cpp
    Distributions/RandomVariable.cpp
)

#CreateCompileGroup(
#    MODELING_SUNDIALS_MODELS
#    "ODE integrators using Sundials"
#    muqModeling
#    "MODELING_CORE"
#    "EIGEN3;BOOST;SUNDIALS"
#    ""
#    ODEData.cpp
#    ODE.cpp
#    RootfindingIVP.cpp
#)
