# define the source files in this directory that we want to include in the
#   generic interface library
SET(SOURCES
    WorkPiece.cpp
  )

##########################################
#  Objects to build

# Define the library
ADD_LIBRARY(muqModelling ${SOURCES})
TARGET_LINK_LIBRARIES(muqModelling ${MUQ_LINK_LIBS})

install(TARGETS muqModelling
    EXPORT MUQDepends
    LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
