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
    SAMPLING_ALGORITHM
    "Implement sampling algorithms such as Monte Carlo, importance sampling, and Markov chain Monte Carlo"
    muqSamplingAlgorithms
    "MODELING_CORE"
    ""
    ""
    SamplingAlgorithm.cpp
    MonteCarlo.cpp
    ImportanceSampling.cpp
    MCMC.cpp

    TransitionKernel.cpp
    MCKernel.cpp
    IPKernel.cpp
    MHKernel.cpp

    SamplingState.cpp

    SamplingProblem.cpp
)
