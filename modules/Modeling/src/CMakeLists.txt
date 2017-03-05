# define the source files in this directory that we want to include in the
#   generic interface library
SET(SOURCES
    Core/WorkPiece.cpp

    Core/WorkGraph.cpp
    Core/WorkGraphNode.cpp
    Core/WorkGraphEdge.cpp
  )

##########################################
#  Objects to build

# Define the library
ADD_LIBRARY(muqModelingCore ${SOURCES})
TARGET_LINK_LIBRARIES(muqModelingCore ${MUQ_LINK_LIBS})

install(TARGETS muqModelingCore
    EXPORT MUQDepends
    LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
