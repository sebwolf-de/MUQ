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

/*
  AdmissableNeighbor.ValidNeighbor
  --------------------------------

  Purpose:
  Make sure the MultiIndex class will return valid admissable neighbors.

  Test:
  Create MultiIndexSet with dimension of 2 and max order of 1. Add an admissable
  neighbor to index (2,0), and check to make sure that index is an admissable
  neighbor to the MultiIndexSet.

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

/*
  AdmissableNeighbor.UndefinedNeighbor
  ------------------------------------

  Purpose:
  Make sure that MultiIndex class will not return a valid neighbor for indices
  that do not have defined admissable neighbors.

  Test:
  Create MultiIndexSet with dimension of 2 and max order of 1. Add an admissable
  neighbor to index (2,1) and, check to make sure that index is not returned as
  an admissable neighbor to the MultiIndexSet - missing a required neighboring
  index (2,0).

  3 |                          3 |
    |                            |
  2 |                          2 |
    |                  --->      |
  1 | x   x   o                1 | x   x   o
    |                            |
  0 | x   x                    0 | x   x
    -----------------            -----------------
      0   1   2   3                0   1   2   3
*/
TEST(AdmissableNeighbor, UndefinedNeighbor)
{
  // MultiIndexSet - the "square".
  shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1);

  // Create MultiIndex for testing against the square.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,2);
  multi.SetValue(1,1);

  // Check the result of IsAdmissable().
  bool admiss = indexFamily->IsAdmissible(multi);
  EXPECT_FALSE( admiss);
}

/*
  AdmissableNeighbor.OutsideMaxOrder
  ----------------------------------

  Purpose:
  Make sure that MultiIndex class will not return a valid neighbor for indices
  that are outside the max order limit.

  Test:
  Create MultiIndexSet with dimension of 2 and order of 3, and create a max
  index limiter of 3. Add an admissable neighbor to index (4,0), and check to
  make sure that index is not returned as an admissable neighbor to the
  MultiIndexSet.

  4 |                           4 |
    -----------------             -----------------
  3 | x   x   x   x |           3 | x   x   x   x |
    |               |             |               |
  2 | x   x   x   x |           2 | x   x   x   x |
    |               |    --->     |               |
  1 | x   x   x   x |           1 | x   x   x   x |
    |               |             |               |
  0 | x   x   x   x | o         0 | x   x   x   x | o
    --------------------          --------------------
      0   1   2   3   4             0   1   2   3   4

*/
TEST(AdmissableNeighbor, OutsideMaxOrder)
{
  // Max order limit of 3.
  shared_ptr<MultiIndexLimiter> limiter = MaxOrderLimiter(3);

  // MultiIndexSet - the "square".
  shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 3, limiter);

  // Create MultiIndex for testing against the square.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,4);

  // Check the result of IsAdmissable().
  bool admiss = indexFamily->IsAdmissible(multi);
  EXPECT_FALSE(admiss);
}

/*
  AdmissableNeighbor.OutsideTotalOrder
  ------------------------------------

  Purpose:
  Make sure that MultiIndex class will not return a valid neighbor for indices
  that are outside the total order.

  Test:
  Create MultiIndexSet with dimension of 2 and max order of 3, and create a
  total index limiter of 3. Add two admissable neighbors to indices (2,2) and
  (4,0), and check to make sure that both indices are not returned as
  admissable neighbors to the MultiIndexSet.

  4 |                           4 |
    |\                            |\
  3 | x\                        3 | x\
    |    \                        |    \
  2 | x   x\  o                 2 | x   x\  o
    |        \          --->      |        \
  1 | x   x   x\                1 | x   x   x\
    |            \                |            \
  0 | x   x   x   x\  o         0 | x   x   x   x\  o
    --------------------          --------------------
      0   1   2   3   4             0   1   2   3   4
*/
TEST(AdmissableNeighbor, OutsideTotalOrder)
{
  // Total order limit of 3.
  shared_ptr<MultiIndexLimiter> limiter = TotalOrderLimiter(3);

  // MultiIndexSet - the "square".
  shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 3, limiter);

  // Create MultiIndex for testing against the square.
  shared_ptr<MultiIndex> multi1 = make_shared<MultiIndex>(2);
  multi1.SetValue(0,4);

  shared_ptr<MultiIndex> multi2 = make_shared<MultiIndex>(2);
  multi2.SetValue(0,2);
  multi2.SetValue(1,2);

  // Check the result of IsAdmissable().
  bool admiss1 = indexFamily->IsAdmissible(multi1);
  EXPECT_FALSE(admiss1);

  bool admiss2 = indexFamily->IsAdmissible(multi2);
  EXPECT_FALSE(admiss2);
}

/*
  AdmissableNeighbor.AddAdmissibleNeighbor
  ----------------------------------------

  Purpose:
  Make sure that MultiIndex class will return a valid neighbor for admissible
  neighbors that have been added to the MultiIndexSet.

  Test:
  Create MultiIndexSet with dimension of 2 and max order of 1, add an admissible
  neighbor at (2,0), then add an admissable neighbor to index (2,1). Then check
  to make sure that new index is returned as an admissable neighbor to the
  MultiIndexSet.

  3 |                          3 |
    |                            |
  2 |                          2 |
    |                  --->      |
  1 | x   x   o                1 | x   x   x
    |                            |
  0 | x   x   x                0 | x   x   x
    ----------------             ----------------
      0   1   2   3                0   1   2   3
*/
TEST(AdmissableNeighbor, AddAdmissibleNeighbor)
{
  // MultiIndexSet - the "square".
  shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1, limiter);

  // Add forward admissible neighbor to index (1,0).
  shared_ptr<MultiIndexSet> localIndex = (2,0)
  indexFamily->AddForwardNeighbors(localIndex)

  // Create MultiIndex for testing against the square.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,2);
  multi.SetValue(1,1);

  // Check the result of IsAdmissable().
  bool admiss = indexFamily->IsAdmissible(multi);
  EXPECT_TRUE(admiss);
}

/*
  AdmissableNeighbor.ForciblyExpandAdmissibleNeighbors
  ----------------------------------------------------

  Purpose:
  Make sure that MultiIndex class will return a valid neighbor for admissible
  neighbors that have been added to the MultiIndexSet.

  Test:
  Create MultiIndexSet with dimension of 2 and max order of 1, add admissible
  neighbors at (1,1) using the ForciblyExpand function. Add an admissable
  neighbor to index (3,0), and check to make sure that index is returned as an
  admissable neighbor to the MultiIndexSet.

  3 |                         3 |
    |                           |
  2 | x   x   x               2 | x   x   x
    |                  --->     |
  1 | x   x   x               1 | x   x   x
    |                           |
  0 | x   x   x   o           0 | x   x   x   x
    -----------------           -----------------
      0   1   2   3               0   1   2   3
*/
TEST(AdmissableNeighbor, ForciblyExpandAdmissibleNeighbors)
{
  // MultiIndexSet - the "square".
  shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1);

  // Add forward admissible neighbor to index (1,1) using ForciblyExpand.
  shared_ptr<MultiIndexSet> localIndex = (1,1)
  indexFamily->ForciblyExpand(localIndex)

  // Create MultiIndex for testing against the square.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,3);

  // Check the result of IsAdmissable().
  bool admiss = indexFamily->IsAdmissible(multi);
  EXPECT_TRUE(admiss);
}


/*
  ------------------------------------------------------------------------------
                           MultiIndexLimiter Tests:
  ------------------------------------------------------------------------------

  Tests whether different MultiIndex limiter classes function correctly.
*/



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



/*
  MultiIndexLimiter.ValidMaxLimiter
  ---------------------------------

  Purpose:
  Make sure the MultiIndexLimiter class returns true value for index within the
  limiter.

  Test:
  Create MultiIndex at (0,0) and max order limiter of 2. Check to make sure the
  MultiIndex is feasible within the limiter.

    -----------                -----------
  2 |         |              2 |         |
    |         |                |         |
  1 |         |     --->     1 |         |
    |         |                |         |
  0 | x       |              0 | x       |
    -----------                -----------
      0   1   2                  0   1   2
*/
TEST(MultiIndexLimiter, ValidMaxLimiter)
{
  // MultiIndexLimiter.
  shared_ptr<MultiIndexLimiter> limiter = MaxOrderLimiter(2);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);

  // Check the result of IsFeasible().
  bool feasible = limiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);
}

/*
  MultiIndexLimiter.InvalidMaxLimiter
  -----------------------------------

  Purpose:
  Make sure the MultiIndexLimiter class returns false value for index outside of
  the max order limiter.

  Test:
  Create MultiIndex at (3,3) and max order limiter of 2. Check to make sure the
  MultiIndex is not feasible within the limiter.

  3 |           x              3 |           o
    -----------                  -----------
  2 |         |                2 |         |
    |         |                  |         |
  1 |         |       --->     1 |         |
    |         |                  |         |
  0 |         |                0 |         |
    -------------                -------------
    0   1   2   3                0   1   2   3
*/
TEST(MultiIndexLimiter, InvalidMaxLimiter)
{
  // MultiIndexLimiter.
  shared_ptr<MultiIndexLimiter> limiter = MaxOrderLimiter(2);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,3);
  multi.SetValue(1,3);

  // Check the result of IsFeasible().
  bool feasible = limiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);
}
/*

/*
  MultiIndexLimiter.ValidTotalLimiter
  -----------------------------------

  Purpose:
  Make sure the MultiIndexLimiter class returns true value for index within the
  total order limiter.

  Test:
  Create MultiIndex at (0,0) and total order limiter of 2. Check to make sure the
  MultiIndex is feasible within the limiter.


  2 |\                       2 |\
    |  \                       |  \
  1 |    \          --->     1 |    \
    |      \                   |      \
  0 | x      \               0 | x      \
    -----------                -----------
      0   1   2                  0   1   2
*/
TEST(MultiIndexLimiter, ValidTotalLimiter)
{
  // MultiIndexLimiter.
  shared_ptr<MultiIndexLimiter> limiter = TotalOrderLimiter(2);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);

  // Check the result of IsFeasible().
  bool feasible = limiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);
}

/*
  MultiIndexLimiter.InvalidTotalLimiter
  -------------------------------------

  Purpose:
  Make sure the MultiIndexLimiter class returns false value for index outside of
  the total order limiter.

  Test:
  Create MultiIndex at (2,2) and total order limiter of 2. Check to make sure the
  MultiIndex is not feasible within the limiter.


2 |\        x              2 |\        o
  |  \                       |  \
1 |    \          --->     1 |    \
  |      \                   |      \
0 |        \               0 |        \
  -----------                -----------
    0   1   2                  0   1   2
*/
TEST(MultiIndexLimiter, InvalidTotalLimiter)
{
  // MultiIndexLimiter.
  shared_ptr<MultiIndexLimiter> limiter = TotalOrderLimiter(2);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,2);
  multi.SetValue(1,2);

  // Check the result of IsFeasible().
  bool feasible = limiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);
}
/*

/*
  MultiIndexLimiter.ValidAndLimiter
  ---------------------------------

  Purpose:
  Make sure the MultiIndexLimiter class returns true value for index within the
  and limiter.

  Test:
  Create MultiIndex at (0,0) with total order limiter of 2 and max order limiter
  of 2. Check to make sure the MultiIndex is feasible within both limiters.

    ------------                -----------
  2 |\         |              2 |\         |
    |  \       |                |  \       |
  1 |    \     |     --->     1 |    \     |
    |      \   |                |      \   |
  0 | x      \ |              0 | x      \ |
    -----------                 ------------
      0   1   2                  0   1   2
*/
TEST(MultiIndexLimiter, ValidAndLimiter)
{
  // MultiIndexLimiter.
  shared_ptr<MultiIndexLimiter> totalLimiter = TotalOrderLimiter(2);
  shared_ptr<MultiIndexLimiter> MaxLimiter = MaxOrderLimiter(2);

  // Combine into AndLimiter.
  shared_ptr<MultiIndexLimiter> andLimiter = AndLimiter(totalLimiter, maxLimiter)

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);

  // Check the result of IsFeasible().
  bool feasible = andLimiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);
}

/*
  MultiIndexLimiter.InvalidAndLimiter1
  ------------------------------------------

  Purpose:
  Make sure the MultiIndexLimiter class returns false value for index outside of
  the total order limiter.

  Test:
  Create MultiIndex at (2,2) with total order limiter of 2 and max order limiter
  of 2. Check to make sure the MultiIndex is not feasible within both limiters.

  ------------                -----------
2 |\        x|              2 |\        o|
  |  \       |                |  \       |
1 |    \     |     --->     1 |    \     |
  |      \   |                |      \   |
0 |        \ |              0 |        \ |
  -----------                 ------------
    0   1   2                  0   1   2
*/
TEST(MultiIndexLimiter, InvalidAndLimiter1)
{
  // MultiIndexLimiter.
  shared_ptr<MultiIndexLimiter> totalLimiter = TotalOrderLimiter(2);
  shared_ptr<MultiIndexLimiter> MaxLimiter = MaxOrderLimiter(2);

  // Combine into AndLimiter.
  shared_ptr<MultiIndexLimiter> andLimiter = AndLimiter(totalLimiter, maxLimiter)

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,2);
  multi.SetValue(1,2);

  // Check the result of IsFeasible().
  bool feasible = andLimiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);
}
/*

/*
  MultiIndexLimiter.InvalidAndLimiter2
  ------------------------------------------

  Purpose:
  Make sure the MultiIndexLimiter class returns false value for index outside of
  the total order limiter.

  Test:
  Create MultiIndex at (3,2) with total order limiter of 2 and max order limiter
  of 2. Check to make sure the MultiIndex is not feasible within both limiters.

  ------------                -----------
2 |\         | x            2 |\         |o
  |  \       |                |  \       |
1 |    \     |     --->     1 |    \     |
  |      \   |                |      \   |
0 |        \ |              0 |        \ |
  -----------                 ------------
    0   1   2                  0   1   2
*/
TEST(MultiIndexLimiter, InvalidAndLimiter2)
{
  // MultiIndexLimiter.
  shared_ptr<MultiIndexLimiter> totalLimiter = TotalOrderLimiter(2);
  shared_ptr<MultiIndexLimiter> MaxLimiter = MaxOrderLimiter(2);

  // Combine into AndLimiter.
  shared_ptr<MultiIndexLimiter> andLimiter = AndLimiter(totalLimiter, maxLimiter)

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,3);
  multi.SetValue(1,2);

  // Check the result of IsFeasible().
  bool feasible = andLimiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);
}
/*

/*
  MultiIndexLimiter.ValidOrLimiter1
  ---------------------------------

  Purpose:
  Make sure the OrLimiter class returns True for case where both constituent
  limiters are True.

  Test:
  Create a TotalOrderLimiter of 2 and MaxOrderLimiter of 2, then create an
  OrLimiter with both. Create a MultiIndex at (1,1) and make sure it is
  feasible - should be True for both constituent limiters.
*/
TEST(MultiIndexLimiter, ValidOrLimiter)
{
  // Max and otal limiters.
  shared_ptr<MultiIndexLimiter> totalLimiter = TotalOrderLimiter(2);
  shared_ptr<MultiIndexLimiter> maxLimiter = MaxOrderLimiter(2);

  // Combine into OrLimiter.
  shared_ptr<MultiIndexLimiter> orLimiter = OrLimiter(totalLimiter, maxLimiter)

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,1);
  multi.SetValue(1,1);

  // Check the result of IsFeasible().
  bool feasible = orLimiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);
}
/*

/*
  MultiIndexLimiter.ValidOrLimiter2
  ---------------------------------

  Purpose:
  Make sure the OrLimiter class returns True for case where one constituent
  limiter is True, and the other is False.

  Test:
  Create a TotalOrderLimiter of 2 and MaxOrderLimiter of 2, then create an
  OrLimiter with both. Create a MultiIndex at (2,1) and make sure it is
  feasible - should be True for MaxOrderLimiter, but False for TotalOrderLimiter.
*/
TEST(MultiIndexLimiter, ValidOrLimiter2)
{
  // Max and otal limiters.
  shared_ptr<MultiIndexLimiter> totalLimiter = TotalOrderLimiter(2);
  shared_ptr<MultiIndexLimiter> maxLimiter = MaxOrderLimiter(2);

  // Combine into OrLimiter.
  shared_ptr<MultiIndexLimiter> orLimiter = OrLimiter(totalLimiter, maxLimiter)

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,2);
  multi.SetValue(1,1);

  // Check the result of IsFeasible().
  bool feasible = orLimiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);
}
/*

/*
  MultiIndexLimiter.InvalidOrLimiter
  ----------------------------------

  Purpose:
  Make sure the OrLimiter class returns False for case where both constituent
  limiters are False.

  Test:
  Create a TotalOrderLimiter of 2 and MaxOrderLimiter of 2, then create an
  OrLimiter with both. Create a MultiIndex at (2,3) and make sure it is NOT
  feasible - should be False for both constituent limiters.
*/
TEST(MultiIndexLimiter, InvalidOrLimiter)
{
  // Max and otal limiters.
  shared_ptr<MultiIndexLimiter> totalLimiter = TotalOrderLimiter(2);
  shared_ptr<MultiIndexLimiter> maxLimiter = MaxOrderLimiter(2);

  // Combine into OrLimiter.
  shared_ptr<MultiIndexLimiter> orLimiter = OrLimiter(totalLimiter, maxLimiter)

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,2);
  multi.SetValue(1,3);

  // Check the result of IsFeasible().
  bool feasible = orLimiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);
}
/*

/*
  MultiIndexLimiter.ValidXorLimiter
  ---------------------------------

  Purpose:
  Make sure the XorLimiter class returns False for case where one constituent
  limiter is True, and the other is False.

  Test:
  Create a TotalOrderLimiter of 2 and MaxOrderLimiter of 2, then create an
  XorLimiter with both. Create a MultiIndex at (2,1) and make sure it is
  feasible - should be True for MaxOrderLimiter, False for TotalOrderLimiter ->
  XOR = True.
*/
TEST(MultiIndexLimiter, ValidXorLimiter)
{
  // Max and otal limiters.
  shared_ptr<MultiIndexLimiter> totalLimiter = TotalOrderLimiter(2);
  shared_ptr<MultiIndexLimiter> maxLimiter = MaxOrderLimiter(2);

  // Combine into XorLimiter.
  shared_ptr<MultiIndexLimiter> XorLimiter = XorLimiter(totalLimiter, maxLimiter)

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,2);
  multi.SetValue(1,1);

  // Check the result of IsFeasible().
  bool feasible = XorLimiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);
}
/*

/*
  MultiIndexLimiter.InvalidXorLimiter1
  ------------------------------------

  Purpose:
  Make sure the XorLimiter class returns False for case where both constituent
  limiters are True.

  Test:
  Create a TotalOrderLimiter of 2 and MaxOrderLimiter of 2, then create an
  XorLimiter with both. Create a MultiIndex at (1,1) and make sure it is NOT
  feasible - should be True for both constituent limiters -> XOR = False.
*/
TEST(MultiIndexLimiter, InvalidXorLimiter1)
{
  // Max and otal limiters.
  shared_ptr<MultiIndexLimiter> totalLimiter = TotalOrderLimiter(2);
  shared_ptr<MultiIndexLimiter> maxLimiter = MaxOrderLimiter(2);

  // Combine into XorLimiter.
  shared_ptr<MultiIndexLimiter> XorLimiter = XorLimiter(totalLimiter, maxLimiter)

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,1);
  multi.SetValue(1,1);

  // Check the result of IsFeasible().
  bool feasible = XorLimiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);
}
/*

/*
  MultiIndexLimiter.InvalidXorLimiter2
  ------------------------------------

  Purpose:
  Make sure the XorLimiter class returns False for case where both constituent
  limiters are False.

  Test:
  Create a TotalOrderLimiter of 2 and MaxOrderLimiter of 2, then create an
  XorLimiter with both. Create a MultiIndex at (2,3) and make sure it is NOT
  feasible - should be False for both constituent limiters -> XOR = False.
*/
TEST(MultiIndexLimiter, InvalidXorLimiter2)
{
  // Max and otal limiters.
  shared_ptr<MultiIndexLimiter> totalLimiter = TotalOrderLimiter(2);
  shared_ptr<MultiIndexLimiter> maxLimiter = MaxOrderLimiter(2);

  // Combine into XorLimiter.
  shared_ptr<MultiIndexLimiter> XorLimiter = XorLimiter(totalLimiter, maxLimiter)

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,2);
  multi.SetValue(1,3);

  // Check the result of IsFeasible().
  bool feasible = XorLimiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);
}
/*
  MultiIndexLimiter.ValidGeneralLimiter
  -------------------------------------

  Purpose:
  Make sure the GeneralLimiter class returns True value for index within the
  GeneralLimiter.

  Test:
  Create MultiIndexSet with max order 1 and dimension 2, and input it to a
  GeneralLimiter. Create MultiIndex at (0,0) and max order of 2. Check to make
  sure the MultiIndex is feasible within the GeneralLimiter.

  2 |                       2 |
    |                         |
  1 | x   x         --->    1 |
    |                         |
  0 | x   x                 0 | x
    -----------               -----------
      0   1   2                 0   1   2
*/
TEST(MultiIndexLimiter, ValidGeneralLimiter)
{
  // MultiIndexSet - the "square".
  shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1);

  // GeneralLimiter.
  shared_ptr<MultiIndexLimiter> limiter = GeneralLimiter(indexFamily);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);

  // Check the result of IsFeasible().
  bool feasible = limiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);
}

/*
  MultiIndexLimiter.InvalidGeneralLimiter
  ---------------------------------------

  Purpose:
  Make sure the GeneralLimiter class returns False value for index outside the
  GeneralLimiter.

  Test:
  Create MultiIndexSet with max order 1 and dimension 2, and input it to a
  GeneralLimiter. Create MultiIndex at (2,1) and max order of 2. Check to make
  sure the MultiIndex is NOT feasible within the GeneralLimiter.

  2 |                       2 |
    |                         |
  1 | x   x         --->    1 |         o
    |                         |
  0 | x   x                 0 |
    -----------               -----------
      0   1   2                 0   1   2
*/
TEST(MultiIndexLimiter, InvalidGeneralLimiter)
{
  // MultiIndexSet - the "square".
  shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1);

  // GeneralLimiter.
  shared_ptr<MultiIndexLimiter> limiter = GeneralLimiter(indexFamily);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(2);
  multi.SetValue(0,2);
  multi.SetValue(1,1);

  // Check the result of IsFeasible().
  bool feasible = limiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);
}

/*
  MultiIndexLimiter.ValidDimensionLimiter
  ---------------------------------------

  Purpose:
  Make sure the MultiIndexLimiter class returns True value for index within the
  DimensionLimiter.

  Test:
  Create DimensionLimiter with lower dimension of 1 and length of 1, and create
  a MultiIndex of [1 0 0 1]. Check to make sure the MultiIndex is feasible
  within the DimensionLimiter.

  DimensionLimiter(1,1) of [a0  a1  a2  a3] =

*/
TEST(MultiIndexLimiter, ValidDimensionLimiter)
{
  // DimensionLimiter.
  shared_ptr<MultiIndexLimiter> limiter = DimensionLimiter(1,1);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(4);
  multi.SetValue(0,1);
  multi.SetValue(3,1);

  // Check the result of IsFeasible().
  bool feasible = limiter->IsFeasible(multi);
  EXPECT_TRUE(feasible);
}

/*
  MultiIndexLimiter.InvalidDimensionLimiter
  -----------------------------------------

  Purpose:
  Make sure the MultiIndexLimiter class returns False value for index within the
  DimensionLimiter.

  Test:
  Create DimensionLimiter with lower dimension of 1 and length of 1, and create
  a MultiIndex of [0 1 0 1]. Check to make sure the MultiIndex is NOT feasible
  within the DimensionLimiter.

*/
TEST(MultiIndexLimiter, InvalidDimensionLimiter)
{
  // DimensionLimiter.
  shared_ptr<MultiIndexLimiter> limiter = DimensionLimiter(1,1);

  // MultiIndex for testing.
  shared_ptr<MultiIndex> multi = make_shared<MultiIndex>(4);
  multi.SetValue(1,1);
  multi.SetValue(3,1);

  // Check the result of IsFeasible().
  bool feasible = limiter->IsFeasible(multi);
  EXPECT_FALSE(feasible);
}
