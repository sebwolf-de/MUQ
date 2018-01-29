
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
  Make sure that MultiIndex class will return valid admissable neighbors.

  Test:
  Create MultiIndexSet with dimension of 2 and max order of 1, and create an
  index limiter of 3. Add two admissable neighbors to indices (2,0) and (0,2),
  and check that those two indices are admissable neighbors to the
  MultiIndexSet.

*/
TEST(AdmissableNeighbor, ValidNeighbor)
{
  // Create max order limit of 3.
  shared_ptr<MultiIndexLimiter> limiter = MaxOrderLimiter(3)

  // MultiIndexSet - the "square".
  shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1, limiter);

  // MultiIndex
  MultiIndex multi(2);
  multi.SetValue(2,0);
  multi.SetValue(0,2);
  Eigen::RowVectorXu temp1 = multi.GetMulti();

  // Check the result of IsAdmissable().
  bool trueBool = true;
  bool admiss = indexFamily->IsAdmissible(temp1);
  EXPECT_EQ(trueBool, admiss);

}

/*

  AdmissableNeighbor.UndefinedNeighbor
  ------------------------------------

  Purpose:
  Make sure that MultiIndex class will not return a valid neighbor for indices
  that do not have defined admissable neighbors.

  Test:
  Create MultiIndexSet with dimension of 2 and max order of 1, and create an
  index limiter of 3. Add two admissable neighbors to indices (2,1) and (1,2),
  and check that those two indices are not identified admissable neighbors to the
  MultiIndexSet.

*/
TEST(AdmissableNeighbor, UndefinedNeighbor)
{
  // Create max order limit of 3.
  shared_ptr<MultiIndexLimiter> limiter = MaxOrderLimiter(3)

  // MultiIndexSet - the "square".
  shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 1, limiter);

  // Create MultiIndex for testing against the square.
  MultiIndex multi(2);
  multi.SetValue(2,1);
  multi.SetValue(1,2);
  Eigen::RowVectorXu temp2 = multi.GetMulti();

  // Check the result of IsAdmissable().
  bool trueBool = true;
  bool admiss = indexFamily->IsAdmissible(temp2);
  EXPECT_EQ();
}

/*

  AdmissableNeighbor.OutsideMaxOrder
  ----------------------------------

  Purpose:
  Make sure that MultiIndex class will not return a valid neighbor for indices
  that are outside the max order.

  Test:
  Create MultiIndexSet with dimension of 2 and max order of 3, and create an
  index limiter of 3. Add two admissable neighbors to indices (4,0) and (0,4),
  and check that those two indices are not identified admissable neighbors to the
  MultiIndexSet.

*/
TEST(AdmissableNeighbor, OutsideMaxOrder)
{
  // Create max order limit of 3.
  shared_ptr<MultiIndexLimiter> limiter = MaxOrderLimiter(3)

  // MultiIndexSet - the "square".
  shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 3, limiter);

  // Create MultiIndex for testing against the square.
  MultiIndex multi(2);
  multi.SetValue(4,0);
  multi.SetValue(0,4);
  Eigen::RowVectorXu temp3 = multi.GetMulti();

  // Check the result of IsAdmissable().
  bool trueBool = true;
  bool admiss = indexFamily->IsAdmissible(temp3);


  // MaxOrderLimiter?


  EXPECT_EQ();
}

/*

  AdmissableNeighbor.OutsideTotalOrder
  ------------------------------------

  Purpose:
  Make sure that MultiIndex class will not return a valid neighbor for indices
  that are outside the total order.

  Test:
  Create MultiIndexSet with dimension of 2 and max order of 3, and create an
  index limiter of 3. Add three admissable neighbors to indices (2,2), (4,0), and (0,4),
  and check that those indices are not identified as admissable neighbors to the
  MultiIndexSet.

*/
TEST(AdmissableNeighbor, OutsideTotalOrder)
{
  // Create total order limit of 3.
  shared_ptr<MultiIndexLimiter> limiter = TotalOrderLimiter(3)

  // MultiIndexSet - the "square".
  shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(2, 3, limiter);

  // Create MultiIndex for testing against the square.
  MultiIndex multi(3);
  multi.SetValue(4,0);
  multi.SetValue(0,4);
  multi.SetValue(2,2);
  Eigen::RowVectorXu temp3 = multi.GetMulti();

  // Check the result of IsAdmissable().
  bool trueBool = true;
  bool admiss = indexFamily->IsAdmissible(temp3);

  EXPECT_EQ();
}
