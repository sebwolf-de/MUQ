# USAGE:
# CreateCompileGroup(
#     <group_name>
#     <brief description>
#     <other required groups>
#     <library_name>
#     <required dependencies>
#     <optional dependencies>
#     <source file 1>
#     <source file 2>
#     ...
#     <source file N>
# )

CreateCompileGroup(
    MODELING_CORE
    muqModelingCore
    "EIGEN3;BOOST"
    ""
    WorkPiece.cpp
)
