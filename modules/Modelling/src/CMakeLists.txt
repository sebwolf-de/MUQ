# USAGE:
# CreateCompileGroup(
#     <group_name>
#     <library_name> 
#     <requirements>
#     <source file 1>
#     <source file 2>
#     ...
#     <source file N>
# )

CreateCompileGroup(
    MODELING_CORE
    muqModelingCore
    "EIGEN"
    WorkPiece.cpp
)
