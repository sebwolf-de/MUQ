#include <iostream>

#include "gtest/gtest.h"

#include "MUQ/Utilities/MultiIndices/MultiIndex.h"
//#include "MUQ/Utilities/MultiIndices/MultiIndexSet.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexLimiter.h"

using namespace std;

using namespace muq::Utilities;

/*
  MultiIndex Diagrams
  -------------------

    x=member of the MultiIndexSet,
    o=neighbor being tested -> changed to x if neighbor is admissible.
*/

/*
  ------------------------------------------------------------------------------
                            AdmissableNeighbor Tests:
  ------------------------------------------------------------------------------

  Tests whether admissable neighbors are identified correctly for different
  scenarios.
*/

/*
  MultiIndexLimiter.ValidMultiIndexLimiter
  ----------------------------------------

  Purpose:
  Make sure the MultiIndexLimiter class returns true value for index within the
  limiter.

  Test:
  Create MultiIndex at (0,0) and max ordrer limiter of 2. Check to make sure the
  MultiIndex is feasible within the limiter.

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
// TEST(AdmissableNeighbor, ValidNeighbor)
// {
//   // MultiIndexSet - the "square".
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1);

//   // MultiIndex for testing.
//   shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
//   multi.SetValue(0,2);

//   // Check the result of IsAdmissable().
//   bool admiss = indexFamily->IsAdmissible(multi);
//   EXPECT_TRUE(admiss);
// }

// /*
//   AdmissableNeighbor.ValidNeighbor
//   --------------------------------

//   Purpose:
//   Make sure the MultiIndex class will return valid admissable neighbors.

//   Test:
//   Create MultiIndexSet with dimension of 2 and max order of 1. Add an admissable
//   neighbor to index (2,0), and check to make sure that index is an admissable
//   neighbor to the MultiIndexSet.

//   3 |                         3 |
//     |                           |
//   2 |                         2 |
//     |                 --->      |
//   1 | x   x                   1 | x   x
//     |                           |
//   0 | x   x   o               0 | x   x   x
//     ----------------            ----------------
//       0   1   2   3               0   1   2   3
// */
// TEST(AdmissableNeighbor, ValidNeighbor)
// {
//   // MultiIndexSet - the "square".
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1);

//   // MultiIndex for testing.
//   shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
//   multi.SetValue(0,2);

//   // Check the result of IsAdmissable().
//   bool admiss = indexFamily->IsAdmissible(multi);
//   EXPECT_TRUE(admiss);
// }

// /*
//   AdmissableNeighbor.UndefinedNeighbor
//   ------------------------------------

//   Purpose:
//   Make sure that MultiIndex class will not return a valid neighbor for indices
//   that do not have defined admissable neighbors.

//   Test:
//   Create MultiIndexSet with dimension of 2 and max order of 1. Add an admissable
//   neighbor to index (2,1) and, check to make sure that index is not returned as
//   an admissable neighbor to the MultiIndexSet - missing a required neighboring
//   index (2,0).

//   3 |                          3 |
//     |                            |
//   2 |                          2 |
//     |                  --->      |
//   1 | x   x   o                1 | x   x   o
//     |                            |
//   0 | x   x                    0 | x   x
//     -----------------            -----------------
//       0   1   2   3                0   1   2   3
// */
// TEST(AdmissableNeighbor, UndefinedNeighbor)
// {
//   // MultiIndexSet - the "square".
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1);

//   // Create MultiIndex for testing against the square.
//   shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
//   multi.SetValue(0,2);
//   multi.SetValue(1,1);

//   // Check the result of IsAdmissable().
//   bool admiss = indexFamily->IsAdmissible(multi);
//   EXPECT_FALSE( admiss);
// }

// /*
//   AdmissableNeighbor.OutsideMaxOrder
//   ----------------------------------

//   Purpose:
//   Make sure that MultiIndex class will not return a valid neighbor for indices
//   that are outside the max order limit.

//   Test:
//   Create MultiIndexSet with dimension of 2 and order of 3, and create a max
//   index limiter of 3. Add an admissable neighbor to index (4,0), and check to
//   make sure that index is not returned as an admissable neighbor to the
//   MultiIndexSet.

//   4 |                           4 |
//     -----------------             -----------------
//   3 | x   x   x   x |           3 | x   x   x   x |
//     |               |             |               |
//   2 | x   x   x   x |           2 | x   x   x   x |
//     |               |    --->     |               |
//   1 | x   x   x   x |           1 | x   x   x   x |
//     |               |             |               |
//   0 | x   x   x   x | o         0 | x   x   x   x | o
//     --------------------          --------------------
//       0   1   2   3   4             0   1   2   3   4

// */
// TEST(AdmissableNeighbor, OutsideMaxOrder)
// {
//   // Max order limit of 3.
//   shared_ptr<MultiIndexLimiter> limiter = MaxOrderLimiter(3);

//   // MultiIndexSet - the "square".
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 3, limiter);

//   // Create MultiIndex for testing against the square.
//   shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
//   multi.SetValue(0,4);

//   // Check the result of IsAdmissable().
//   bool admiss = indexFamily->IsAdmissible(multi);
//   EXPECT_FALSE(admiss);
// }

// /*
//   AdmissableNeighbor.OutsideTotalOrder
//   ------------------------------------

//   Purpose:
//   Make sure that MultiIndex class will not return a valid neighbor for indices
//   that are outside the total order.

//   Test:
//   Create MultiIndexSet with dimension of 2 and max order of 3, and create a
//   total index limiter of 3. Add two admissable neighbors to indices (2,2) and
//   (4,0), and check to make sure that both indices are not returned as
//   admissable neighbors to the MultiIndexSet.

//   4 |                           4 |
//     |\                            |\
//   3 | x\                        3 | x\
//     |    \                        |    \
//   2 | x   x\  o                 2 | x   x\  o
//     |        \          --->      |        \
//   1 | x   x   x\                1 | x   x   x\
//     |            \                |            \
//   0 | x   x   x   x\  o         0 | x   x   x   x\  o
//     --------------------          --------------------
//       0   1   2   3   4             0   1   2   3   4
// */
// TEST(AdmissableNeighbor, OutsideTotalOrder)
// {
//   // Total order limit of 3.
//   shared_ptr<MultiIndexLimiter> limiter = TotalOrderLimiter(3);

//   // MultiIndexSet - the "square".
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 3, limiter);

//   // Create MultiIndex for testing against the square.
//   shared_ptr<MultiIndex> multi1 = make_shared<MultiIndex>(2);
//   multi1.SetValue(0,4);

//   shared_ptr<MultiIndex> multi2 = make_shared<MultiIndex>(2);
//   multi2.SetValue(0,2);
//   multi2.SetValue(1,2);

//   // Check the result of IsAdmissable().
//   bool admiss1 = indexFamily->IsAdmissible(multi1);
//   EXPECT_FALSE(admiss1);

//   bool admiss2 = indexFamily->IsAdmissible(multi2);
//   EXPECT_FALSE(admiss2);
// }

// /*
//   AdmissableNeighbor.AddAdmissibleNeighbor
//   ----------------------------------------

//   Purpose:
//   Make sure that MultiIndex class will return a valid neighbor for admissible
//   neighbors that have been added to the MultiIndexSet.

//   Test:
//   Create MultiIndexSet with dimension of 2 and max order of 1, add an admissible
//   neighbor at (2,0), then add an admissable neighbor to index (2,1). Then check
//   to make sure that new index is returned as an admissable neighbor to the
//   MultiIndexSet.

//   3 |                          3 |
//     |                            |
//   2 |                          2 |
//     |                  --->      |
//   1 | x   x   o                1 | x   x   x
//     |                            |
//   0 | x   x   x                0 | x   x   x
//     ----------------             ----------------
//       0   1   2   3                0   1   2   3
// */
// TEST(AdmissableNeighbor, AddAdmissibleNeighbor)
// {
//   // MultiIndexSet - the "square".
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1, limiter);

//   // Add forward admissible neighbor to index (1,0).
//   shared_ptr<MultiIndexSet> localIndex = (2,0)
//   indexFamily->AddForwardNeighbors(localIndex)

//   // Create MultiIndex for testing against the square.
//   shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
//   multi.SetValue(0,2);
//   multi.SetValue(1,1);

//   // Check the result of IsAdmissable().
//   bool admiss = indexFamily->IsAdmissible(multi);
//   EXPECT_TRUE(admiss);
// }

// /*
//   AdmissableNeighbor.ForciblyExpandAdmissibleNeighbors
//   ----------------------------------------------------

//   Purpose:
//   Make sure that MultiIndex class will return a valid neighbor for admissible
//   neighbors that have been added to the MultiIndexSet.

//   Test:
//   Create MultiIndexSet with dimension of 2 and max order of 1, add admissible
//   neighbors at (1,1) using the ForciblyExpand function. Add an admissable
//   neighbor to index (3,0), and check to make sure that index is returned as an
//   admissable neighbor to the MultiIndexSet.

//   3 |                         3 |
//     |                           |
//   2 | x   x   x               2 | x   x   x
//     |                  --->     |
//   1 | x   x   x               1 | x   x   x
//     |                           |
//   0 | x   x   x   o           0 | x   x   x   x
//     -----------------           -----------------
//       0   1   2   3               0   1   2   3
// */
// TEST(AdmissableNeighbor, ForciblyExpandAdmissibleNeighbors)
// {
//   // MultiIndexSet - the "square".
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1);

//   // Add forward admissible neighbor to index (1,1) using ForciblyExpand.
//   shared_ptr<MultiIndexSet> localIndex = (1,1)
//   indexFamily->ForciblyExpand(localIndex)

//   // Create MultiIndex for testing against the square.
//   shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
//   multi.SetValue(0,3);

//   // Check the result of IsAdmissable().
//   bool admiss = indexFamily->IsAdmissible(multi);
//   EXPECT_TRUE(admiss);
// }


/*
  ------------------------------------------------------------------------------
                           MultiIndexLimiter Tests:
  ------------------------------------------------------------------------------

  Tests whether different MultiIndex limiter classes function correctly.
*/


/*
IS THIS ACTUALLY GOING TO SAVE ANY SPACE?


class LimiterTest : public testing::TestWithParam<tuple<
                                    shared_ptr<MultiIndexLimiter>, // Limiter
                                    int,                           // MultIndex something
                                    int,
                                    string>>
{

public:
  virtual void SetUp()
  {

  }

  virtual void TearDown()
  {

  }

};

INSTANTIATE_TEST_CASE_P( )


TEST_P(MultiIndexLimiter, ValidMaxLimiter)
{
  shared_ptr<MultiIndexLimiter> limiter = MaxOrderLimiter(2);

  EXPECT_TRUE(LimiterTest(limiter, ));

}
*/


