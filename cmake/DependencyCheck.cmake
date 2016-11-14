# If there is no gtest, turn off all tests
if(NOT MUQ_USE_GTEST)
  set(UtilitiesAndModelling_tests OFF)
  message(WARNING "MUQ_USE_GTEST is OFF. Thus, all tests has been turned off.")
endif(NOT MUQ_USE_GTEST)