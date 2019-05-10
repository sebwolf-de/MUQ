CreateTestGroup(
    SAMPLING_ALGORITHM

    SamplingProblemTests.cpp
    SampleCollectionTests.cpp
    MarkovChainTests.cpp
    MIMCMCTests.cpp
    MLMCMCTests.cpp

    GMHKernelTests.cpp

    #MonteCarloTests.cpp
    #ImportanceSamplingTests.cpp
    MCMCTests.cpp
    CrankNicolsonTests.cpp

    ExpensiveSamplingProblemTests.cpp
)

CreateParallelTestGroup(
    SAMPLING_ALGORITHM

    SamplingStateParallelTests.cpp
    
    DistributedCollectionTests.cpp
    
    GMHKernelParallelTests.cpp

    SingleChainMCMCParallelTests.cpp

    ExpensiveSamplingProblemParallelTests.cpp
)
