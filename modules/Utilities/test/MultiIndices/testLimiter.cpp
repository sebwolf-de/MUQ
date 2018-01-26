#include <iostream>

#include "gtest/gtest.h"

#include "MUQ/Utilities/multiIndex/MultiIndex.h"
#include "MUQ/Utilities/multiIndex/MultiIndexSet.h"
#include "MUQ/Utilities/multiIndex/MultiIndexLimiter.h"

using namespace std;

using namespace muq::Utilities;

/*
  AdmissableNeighbor.ValidNeighbor
  --------------------------------

  Purpose:
  Make sure the MultiIndex class will return valid admissable neighbors.

  Test:
  Create MultiIndexSet with dimension of 2 and max order of 1. Add an admissable
  neighbor to index (2,0), and check to make sure that index is an admissable
  neighbor to the MultiIndexSet.

  MultiIndex Diagram: (x=member of the MultiIndexSet, o=neighbor being tested:
                       changed to x if neighbor is admissible)


  3 |                         3 |
    |                           |
  2 |                         2 |
    |                 --->      |
  1 | x   x                   1 | x   x
    |                           |
  0 | x   x   o               0 | x   x   x
    ----------------            ----------------
      0   1   2   3               0   1   2   3
*/
TEST(AdmissableNeighbor, ValidNeighbor)
{
  // MultiIndexSet - the "square".
  shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,2);

  // Check the result of IsAdmissable().
  bool admiss = indexFamily->IsAdmissible(multi);
  EXPECT_TRUE(admiss);

}
