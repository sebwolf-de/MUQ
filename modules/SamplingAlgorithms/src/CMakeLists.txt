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
    "MODELING_CORE;UTILITIES_CORE"
    ""
    "PARCER"
    SamplingAlgorithm.cpp
    #MonteCarlo.cpp
    #ImportanceSampling.cpp
    SingleChainMCMC.cpp

    SamplingState.cpp
    SampleCollection.cpp
    MarkovChain.cpp

    TransitionKernel.cpp
    #MCKernel.cpp
    #ISKernel.cpp
    MHKernel.cpp
    GMHKernel.cpp
    #MCMCKernel.cpp
    #ParallelMHKernel.cpp
    
    MCMCFactory.cpp

    ThinScheduler.cpp

    MCMCProposal.cpp
    MHProposal.cpp
    AMProposal.cpp
    InverseGammaProposal.cpp
    CrankNicolsonProposal.cpp

    SamplingProblem.cpp
)
