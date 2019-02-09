#include <iostream>

#include "gtest/gtest.h"

#include "MUQ/Approximation/Polynomials/PhysicistHermite.h"
#include "MUQ/Approximation/Polynomials/Legendre.h"

#include "MUQ/Approximation/Quadrature/GaussQuadrature.h"
#include "MUQ/Approximation/Quadrature/FullTensorQuadrature.h"
#include "MUQ/Approximation/Quadrature/SmolyakQuadrature.h"

using namespace muq::Approximation;

///Regression test for a 6-dim mixed var type and mixed order full tensor quadrature grid
TEST(Quadrature, FullTensorSquareTest)
{
  auto legendrePoly = std::make_shared<Legendre>();
  auto gaussLegendre = std::make_shared<GaussQuadrature>(legendrePoly);

  FullTensorQuadrature quad({gaussLegendre, gaussLegendre});
  quad.Compute(5);

  gaussLegendre->Compute(5);
  EXPECT_EQ(std::pow(gaussLegendre->Points().cols(),2), quad.Points().cols());

  EXPECT_NEAR(4.0, quad.Weights().sum(),2e-5);

  Eigen::RowVectorXi multi(2);
  multi << 1, 0;
  quad.Compute(multi);
  EXPECT_EQ(2,quad.Points().cols());

  multi << 0, 1;
  quad.Compute(multi);
  EXPECT_EQ(2,quad.Points().cols());
}

TEST(Quadrature, SmolyakSquareTest)
{
  unsigned int dim = 4;
  unsigned int order = 3;

  auto legendrePoly = std::make_shared<Legendre>();
  auto gaussLegendre = std::make_shared<GaussQuadrature>(legendrePoly);

  SmolyakQuadrature quad(dim, gaussLegendre);
  quad.Compute(order);

  FullTensorQuadrature quad2(dim, gaussLegendre);
  quad2.Compute(order);

  EXPECT_GT(quad2.Points().cols(), quad.Points().cols());

  EXPECT_NEAR(std::pow(2.0,dim), quad.Weights().sum(), 1e-12);
}
