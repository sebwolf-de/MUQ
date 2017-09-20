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
    MODELING_CORE_PYTHON
    "Python wrappers for the basic building blocks of the Modeling library."
    "pymuqModeling"
    "MODELING_CORE"
    "PYTHON"
    ""
    WorkPiecePython.cpp
)

CreateCompileGroup(
    MODELING_DOLFIN
    "Wrappers allowing the construction of WorkPieces Fenics wrappers for the basic building blocks of the Modeling library."
    "pymuqModeling"
    "MODELING_CORE_PYTHON"
    "PYTHON;DOLFIN"
    ""
    FenicsPiece.cpp
)
