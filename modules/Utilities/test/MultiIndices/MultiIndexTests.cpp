
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

TEST(UtilitiesMultiIndex, SetValue)
{

  // Start off with [0, 7]
  MultiIndex multi(2);
  multi.SetValue(1,7);

  EXPECT_EQ(0, multi.GetValue(0));
  EXPECT_EQ(7, multi.GetValue(1));

  EXPECT_EQ(7, multi.Max());
  EXPECT_EQ(7, multi.Sum());
  EXPECT_EQ(2, multi.GetLength());

  // [2, 7]
  multi.SetValue(0,2);
  EXPECT_EQ(2, multi.GetValue(0));
  EXPECT_EQ(7, multi.GetValue(1));
  EXPECT_EQ(7, multi.Max());
  EXPECT_EQ(9, multi.Sum());

  // [2,8]
  multi.SetValue(1,8);
  EXPECT_EQ(2, multi.GetValue(0));
  EXPECT_EQ(8, multi.GetValue(1));
  EXPECT_EQ(8, multi.Max());
  EXPECT_EQ(10, multi.Sum());

  // [2,3]
  multi.SetValue(1,3);
  EXPECT_EQ(2, multi.GetValue(0));
  EXPECT_EQ(3, multi.GetValue(1));
  EXPECT_EQ(3, multi.Max());
  EXPECT_EQ(5, multi.Sum());

 // [2, 0]
  multi.SetValue(1,0);
  EXPECT_EQ(2, multi.GetValue(0));
  EXPECT_EQ(0, multi.GetValue(1));
  EXPECT_EQ(2, multi.Max());
  EXPECT_EQ(2, multi.Sum());
}
