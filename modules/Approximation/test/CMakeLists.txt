
message(STATUS "ADDING APPROXIMATION TESTS!")

CreateTestGroup(
  APPROXIMATION_GP
  GaussianProcesses/KernelTests.cpp
)



