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

    LinearSDE.cpp
)

if(MUQ_USE_DOLFIN)
  CreateCompileGroup(     
      MODELING_DOLFIN
      "Wrappers allowing the construction of WorkPieces Fenics wrappers for the basic building blocks of the Modeling library."
      "pymuqModeling"
      "MODELING_CORE_PYTHON"
      "PYTHON;DOLFIN"
      ""
      Dolfin/DolfinPiece.cpp
  )

endif()
