


set(MUQ_DEPENDS "${MUQ_DEPENDS};EIGEN" CACHE INTERNAL "List of Dependencies that MUQ will need.")


set(MUQ_TARGETS "${MUQ_TARGETS};muqModelingCore" CACHE INTERNAL "List of MUQ libraries to build.")
SET(muqModelingCore_SOURCES

    modules/Modelling/src/WorkPiece.cpp

    CACHE INTERNAL "List of source files for the muqModelingCore library."
)
