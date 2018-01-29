
#include <iostream>

#include "gtest/gtest.h"

#include "MUQ/Utilities/MultiIndices/MultiIndex.h"

using namespace std;
using namespace muq::Utilities;


TEST(UtilitiesMultiIndex, DefaultConstructor)
{
  MultiIndex multi(2);

  EXPECT_EQ(2,multi.GetLength());
  EXPECT_EQ(0,multi.GetValue(0));
  EXPECT_EQ(0,multi.GetValue(1));
  EXPECT_EQ(0,multi.Max());
  EXPECT_EQ(0,multi.Sum());
}

TEST(UtilitiesMultiIndex, RowConstructor)
{
  Eigen::RowVectorXi temp(3);
  temp(0) = 1;
  temp(1) = 5;
  temp(2) = 0;

  MultiIndex multi(temp);

  EXPECT_EQ(3,multi.GetLength());
  EXPECT_EQ(1,multi.GetValue(0));
  EXPECT_EQ(5,multi.GetValue(1));
  EXPECT_EQ(5,multi.Max());
  EXPECT_EQ(6,multi.Sum());
}

TEST(UtilitiesMultiIndex, InitializerListConstructor)
{
  MultiIndex multi{1,5,3,0};

  EXPECT_EQ(4,multi.GetLength());
  EXPECT_EQ(1,multi.GetValue(0));
  EXPECT_EQ(5,multi.GetValue(1));
  EXPECT_EQ(3,multi.GetValue(2));
  EXPECT_EQ(0,multi.GetValue(3));

  EXPECT_EQ(5,multi.Max());
  EXPECT_EQ(9,multi.Sum());
}

//
// TEST(UtilitiesMultiIndex, SetValue)
// {
//   MultiIndex multi(2);
//   multi.SetValue(0,1);
//   multi.SetValue(1,2);
//   Eigen::RowVectorXu temp3 = multi.GetMulti();
//
//   EXPECT_EQ(2,multi.GetDimension());
//   EXPECT_EQ(1,temp3(0));
//   EXPECT_EQ(2,temp3(1));
//   EXPECT_EQ(2,multi.GetMax());
// }
//
// TEST(UtilitiesMultiIndexFactory, FullTensor)
// {
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(5, 5);
//
//   Eigen::RowVectorXu index0    = indexFamily->IndexToMulti(0);
//   Eigen::RowVectorXu index28   = indexFamily->IndexToMulti(28);
//   Eigen::RowVectorXu index837  = indexFamily->IndexToMulti(837);
//   Eigen::RowVectorXu indexLast = indexFamily->IndexToMulti(7775);
//
//   Eigen::RowVectorXu index0Correct(5);
//
//   index0Correct << 0, 0, 0, 0, 0;
//   Eigen::RowVectorXu index28Correct(5);
//   index28Correct << 0, 0, 0, 4, 4;
//   Eigen::RowVectorXu index837Correct(5);
//   index837Correct << 0, 3, 5, 1, 3;
//   Eigen::RowVectorXu indexLastCorrect(5);
//   indexLastCorrect << 5, 5, 5, 5, 5;
//
//   EXPECT_TRUE(MatrixEqual(index0, index0Correct));
//   EXPECT_TRUE(MatrixEqual(index28, index28Correct));
//   EXPECT_TRUE(MatrixEqual(index837, index837Correct));
//   EXPECT_TRUE(MatrixEqual(indexLast, indexLastCorrect));
//
// }
//
//
// TEST(UtilitiesMultiIndexFactory, TotalOrder)
// {
//   //Test the number of a larger family
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateTotalOrder(7, 6);
//
//   EXPECT_EQ(1716u, indexFamily->GetNumberOfIndices());
//
//   Eigen::RowVectorXu index0    = indexFamily->IndexToMulti(0);
//   Eigen::RowVectorXu index49   = indexFamily->IndexToMulti(49);
//   Eigen::RowVectorXu index1647 = indexFamily->IndexToMulti(1647);
//   Eigen::RowVectorXu indexLast = indexFamily->IndexToMulti(1715);
//
//   Eigen::RowVectorXu index0Correct(7);
//   index0Correct << 0, 0, 0, 0, 0, 0, 0;
//   Eigen::RowVectorXu index49Correct(7);
//   index49Correct << 0, 0, 0, 0, 2, 0, 0;
//
//   Eigen::RowVectorXu index1647Correct(7);
//   index1647Correct << 3, 0, 2, 0, 0, 0, 1;
//
//   Eigen::RowVectorXu indexLastCorrect(7);
//   indexLastCorrect << 6, 0, 0, 0, 0, 0, 0;
//
//
//   EXPECT_TRUE(MatrixEqual(index0, index0Correct));
//   EXPECT_TRUE(MatrixEqual(index49, index49Correct));
//   EXPECT_TRUE(MatrixEqual(index1647, index1647Correct));
//   EXPECT_TRUE(MatrixEqual(indexLast, indexLastCorrect));
//
// }
//
//
// TEST(UtilitiesMultiIndexFactory, Hyperbolic)
// {
//   //Test the number of a larger family
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateHyperbolic(7, 6, 0.5);
//
//   EXPECT_EQ(106, indexFamily->GetNumberOfIndices());
//
//   Eigen::RowVectorXu index0    = indexFamily->IndexToMulti(0);
//   Eigen::RowVectorXu index49   = indexFamily->IndexToMulti(49);
//   Eigen::RowVectorXu indexLast = indexFamily->IndexToMulti(105);
//
//   Eigen::RowVectorXu index0Correct(7);
//   index0Correct << 0, 0, 0, 0, 0, 0, 0;
//   Eigen::RowVectorXu index49Correct(7);
//   index49Correct << 0, 0, 1, 0, 2, 0, 0;
//   Eigen::RowVectorXu indexLastCorrect(7);
//   indexLastCorrect << 6, 0, 0, 0, 0, 0, 0;
//
//   EXPECT_TRUE(MatrixEqual(index0, index0Correct));
//   EXPECT_TRUE(MatrixEqual(index49, index49Correct));
//   EXPECT_TRUE(MatrixEqual(indexLast, indexLastCorrect));
//
// }
//
//
//
// TEST(UtilitiesMultiIndexFactory, LimitedTotalOrder)
// {
//   //Test the number of a larger family
//   auto limiter = make_shared<TotalOrderLimiter>(3);
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateTotalOrder(2, 5,0,limiter);
//
//   EXPECT_EQ(10u, indexFamily->GetNumberOfIndices());
//
//   Eigen::MatrixXu trueMultis(10,2);
//   trueMultis << 0, 0,
//                 0, 1,
//                 0, 2,
//                 0, 3,
//                 1, 0,
//                 1, 1,
//                 1, 2,
//                 2, 0,
//                 2, 1,
//                 3, 0;
//
//   EXPECT_TRUE(MatrixEqual(trueMultis,indexFamily->GetAllMultiIndices()));
// }
//
//
// ///Test the first, last, and a couple of middle multi-indices from the total order family.
// TEST(UtilitiesMultiIndexFactory, TotalOrderWithMin)
// {
//   //Test the number of a larger family
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateTotalOrder(7,8,4);
//
//   EXPECT_EQ(6315u, indexFamily->GetNumberOfIndices());
//
//   Eigen::RowVectorXu index0    = indexFamily->IndexToMulti(0);
//   Eigen::RowVectorXu index53   = indexFamily->IndexToMulti(53);
//   Eigen::RowVectorXu index5461 = indexFamily->IndexToMulti(5461);
//   Eigen::RowVectorXu indexLast = indexFamily->IndexToMulti(6314);
//
//
//   Eigen::RowVectorXu index0Correct(7);
//   index0Correct << 0, 0, 0, 0, 0, 0, 4;
//
//   Eigen::RowVectorXu index53Correct(7);
//   index53Correct << 0, 0, 0, 0, 1, 3, 3;
//
//   Eigen::RowVectorXu index5461Correct(7);
//   index5461Correct << 2, 3, 0, 1, 0, 0, 1;
//
//   Eigen::RowVectorXu indexLastCorrect(7);
//   indexLastCorrect << 8, 0, 0, 0, 0, 0, 0;
//
//
//   EXPECT_TRUE(MatrixEqual(index0, index0Correct));
//   EXPECT_TRUE(MatrixEqual(index53, index53Correct));
//   EXPECT_TRUE(MatrixEqual(index5461, index5461Correct));
//   EXPECT_TRUE(MatrixEqual(indexLast, indexLastCorrect));
//
// }
//
// ///Test the first, last, and a couple of middle multi-indices from the total order family, non-uniform.
// TEST(UtilitiesMultiIndexFactory, NonUniformTensor)
// {
//   //Test the number of a larger family
//   Eigen::RowVectorXu orders(5);
//
//   orders << 1, 2, 4, 3, 4;
//   shared_ptr<MultiIndexSet> indexFamily = MultiIndexFactory::CreateFullTensor(orders);
//
//   EXPECT_EQ(600u, indexFamily->GetNumberOfIndices());
//
//   Eigen::RowVectorXu index0    = indexFamily->IndexToMulti(0);
//   Eigen::RowVectorXu index146  = indexFamily->IndexToMulti(146);
//   Eigen::RowVectorXu index316  = indexFamily->IndexToMulti(316);
//   Eigen::RowVectorXu indexLast = indexFamily->IndexToMulti(599);
//
//   Eigen::RowVectorXu index0Correct(5);
//   index0Correct << 0, 0, 0, 0, 0;
//   Eigen::RowVectorXu index146Correct(5);
//   index146Correct << 0, 1, 2, 1, 1;
//
//   Eigen::RowVectorXu index316Correct(5);
//   index316Correct << 1, 0, 0, 3, 1;
//
//   Eigen::RowVectorXu indexLastCorrect(5);
//   indexLastCorrect << 1, 2, 4, 3, 4;
//
//
//   EXPECT_TRUE(MatrixEqual(index0, index0Correct));
//   EXPECT_TRUE(MatrixEqual(index146, index146Correct));
//   EXPECT_TRUE(MatrixEqual(index316, index316Correct));
//   EXPECT_TRUE(MatrixEqual(indexLast, indexLastCorrect));
//
// }
