CreateTestGroup(
    SAMPLING_ALGORITHM

    SamplingProblemTests.cpp
    SampleCollectionTests.cpp
    MarkovChainTests.cpp

    GMHKernelTests.cpp

    #MonteCarloTests.cpp
    #ImportanceSamplingTests.cpp
    MCMCTests.cpp
    CrankNicolsonTests.cpp
)

CreateParallelTestGroup(
    SAMPLING_ALGORITHM

    SamplingStateParallelTests.cpp
    
    GMHKernelParallelTests.cpp

    SingleChainMCMCParallelTests.cpp
)
