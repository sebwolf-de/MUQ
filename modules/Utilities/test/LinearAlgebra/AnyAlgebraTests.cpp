#include <gtest/gtest.h>

#include "MUQ/Utilities/LinearAlgebra/AnyAlgebra.h"

using namespace muq::Utilities;

TEST(AnyAlgebraTests, Size) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // doubles should have size 1
  EXPECT_EQ(alg->Size(2.0), 1);

  // check for the Eigen::Vector's
  const Eigen::Vector2d test2(2.0, 3.0);
  EXPECT_EQ(alg->Size(test2), 2);
  const Eigen::Vector3d test3(2.0, 3.0, 4.0);
  EXPECT_EQ(alg->Size(test3), 3);
  const Eigen::Vector4d test4(2.0, 3.0, 4.0, 5.0);
  EXPECT_EQ(alg->Size(test4), 4);
  Eigen::VectorXd test(6);
  test << 2.0, 3.0, 4.0, 5.0, 6.0, 7.0;
  EXPECT_EQ(alg->Size(test), 6);

  // check for Eigen::Matrix's
  const Eigen::Matrix2d mat2 = Eigen::Matrix2d::Constant(2.0);
  EXPECT_EQ(alg->Size(mat2), 4); // total number of entries
  EXPECT_EQ(alg->Size(mat2, 0), 2); // number of rows
  EXPECT_EQ(alg->Size(mat2, 1), 2); // number of cols
  const Eigen::Matrix3d mat3 = Eigen::Matrix3d::Constant(2.0);
  EXPECT_EQ(alg->Size(mat3), 9); // total number of entries
  EXPECT_EQ(alg->Size(mat3, 0), 3); // number of rows
  EXPECT_EQ(alg->Size(mat3, 1), 3); // number of cols
  const Eigen::Matrix4d mat4 = Eigen::Matrix4d::Constant(2.0);
  EXPECT_EQ(alg->Size(mat4), 16); // total number of entries
  EXPECT_EQ(alg->Size(mat4, 0), 4); // number of rows
  EXPECT_EQ(alg->Size(mat4, 1), 4); // number of cols
  const Eigen::MatrixXd mat = Eigen::MatrixXd::Constant(4, 5, 2.0);
  EXPECT_EQ(alg->Size(mat), 20); // total number of entries
  EXPECT_EQ(alg->Size(mat, 0), 4); // number of rows
  EXPECT_EQ(alg->Size(mat, 1), 5); // number of cols

  // if muq was compiled with Sundials, check for Sundials types
#if MUQ_HAS_SUNDIALS==1
  N_Vector vec = N_VNew_Serial(8);
  for( unsigned int i=0; i<8; ++i ) { NV_Ith_S(vec, i) = 3.0*i; }

  EXPECT_EQ(alg->Size(vec), 8);

  N_VDestroy(vec);
#endif
}

TEST(AnyAlgebraTests, AccessElement) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test scalar
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(-4.0, 0)), -4.0);

  // test Eigen::Vector's
  const Eigen::Vector2d test2d(2.0, 4.0);
  const Eigen::Vector2f test2f(-32.0, 12.5);
  const Eigen::Vector2i test2i(-5, 4);
  for( unsigned int i=0; i<2; ++i ) {
    EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(test2d, i)), test2d(i));
    EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(test2f, i)), test2f(i));
    EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(test2i, i)), test2i(i));
  }

  const Eigen::Vector3d test3d(3.0, 4.0, -1.3);
  const Eigen::Vector3f test3f(-33.0, 13.5, 0.1);
  const Eigen::Vector3i test3i(-5, 4, 8);
  for( unsigned int i=0; i<3; ++i ) {
    EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(test3d, i)), test3d(i));
    EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(test3f, i)), test3f(i));
    EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(test3i, i)), test3i(i));
  }

  const Eigen::Vector4d test4d(4.0, 4.0, -1.4, 1.3);
  const Eigen::Vector4f test4f(-44.0, 14.5, 0.1, 3.6);
  const Eigen::Vector4i test4i(-5, 4, 8, 3);
  for( unsigned int i=0; i<4; ++i ) {
    EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(test4d, i)), test4d(i));
    EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(test4f, i)), test4f(i));
    EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(test4i, i)), test4i(i));
  }

  const Eigen::VectorXd testXd = Eigen::VectorXd::Random(13);
  const Eigen::VectorXf testXf = Eigen::VectorXf::Random(13);
  const Eigen::VectorXi testXi = Eigen::VectorXi::Random(13);
  for( unsigned int i=0; i<13; ++i ) {
    EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(testXd, i)), testXd(i));
    EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(testXf, i)), testXf(i));
    EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(testXi, i)), testXi(i));
  }

  // test Eigen::Matrix's
  const Eigen::Matrix2d mat2d = Eigen::Matrix2d::Random();
  const Eigen::Matrix2f mat2f = Eigen::Matrix2f::Random();
  const Eigen::Matrix2i mat2i = Eigen::Matrix2i::Random();
  for( unsigned int i=0; i<2; ++i ) {
    for( unsigned int j=0; j<2; ++j ) {
      EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(mat2d, i, j)), mat2d(i, j));
      EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(mat2f, i, j)), mat2f(i,j ));
      EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(mat2i, i, j)), mat2i(i, j));
    }
  }

  const Eigen::Matrix3d mat3d = Eigen::Matrix3d::Random();
  const Eigen::Matrix3f mat3f = Eigen::Matrix3f::Random();
  const Eigen::Matrix3i mat3i = Eigen::Matrix3i::Random();
  for( unsigned int i=0; i<3; ++i ) {
    for( unsigned int j=0; j<3; ++j ) {
      EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(mat3d, i, j)), mat3d(i, j));
      EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(mat3f, i, j)), mat3f(i,j ));
      EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(mat3i, i, j)), mat3i(i, j));
    }
  }

  const Eigen::Matrix4d mat4d = Eigen::Matrix4d::Random();
  const Eigen::Matrix4f mat4f = Eigen::Matrix4f::Random();
  const Eigen::Matrix4i mat4i = Eigen::Matrix4i::Random();
  for( unsigned int i=0; i<4; ++i ) {
    for( unsigned int j=0; j<4; ++j ) {
      EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(mat4d, i, j)), mat4d(i, j));
      EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(mat4f, i, j)), mat4f(i,j ));
      EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(mat4i, i, j)), mat4i(i, j));
    }
  }

  const Eigen::MatrixXd matd = Eigen::MatrixXd::Random(5, 8);
  const Eigen::MatrixXf matf = Eigen::MatrixXf::Random(5, 8);
  const Eigen::MatrixXi mati = Eigen::MatrixXi::Random(5, 8);
  for( unsigned int i=0; i<5; ++i ) {
    for( unsigned int j=0; j<8; ++j ) {
      EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(matd, i, j)), matd(i, j));
      EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->AccessElement(matf, i, j)), matf(i,j ));
      EXPECT_EQ(boost::any_cast<int const>(alg->AccessElement(mati, i, j)), mati(i, j));
    }
  }

  // if muq was compiled with Sundials, check for Sundials types
#if MUQ_HAS_SUNDIALS==1
  N_Vector vec = N_VNew_Serial(8);
  for( unsigned int i=0; i<8; ++i ) { NV_Ith_S(vec, i) = 3.0*i; }

  for( unsigned int i=0; i<8; ++i ) {
    EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->AccessElement(vec, i)), NV_Ith_S(vec, i));
  }

  N_VDestroy(vec);
#endif
}

TEST(AnyAlgebraTests, Zero) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar zero
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Zero(typeid(double))), 0.0);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Zero(typeid(float))), 0.0);
  EXPECT_EQ(boost::any_cast<int const>(alg->Zero(typeid(int))), 0);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Zero(typeid(unsigned int))), 0);

  // test the Eigen::Vector zero
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector2d const&>(alg->Zero(typeid(Eigen::Vector2d))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector2f const&>(alg->Zero(typeid(Eigen::Vector2f))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector2i const&>(alg->Zero(typeid(Eigen::Vector2i))).norm(), 0.0);

  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector3d const&>(alg->Zero(typeid(Eigen::Vector3d))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector3f const&>(alg->Zero(typeid(Eigen::Vector3f))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector3i const&>(alg->Zero(typeid(Eigen::Vector3i))).norm(), 0.0);

  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector4d const&>(alg->Zero(typeid(Eigen::Vector4d))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector4f const&>(alg->Zero(typeid(Eigen::Vector4f))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Vector4i const&>(alg->Zero(typeid(Eigen::Vector4i))).norm(), 0.0);

  const Eigen::VectorXd vecd = boost::any_cast<Eigen::VectorXd const&>(alg->Zero(typeid(Eigen::VectorXd), 13));
  const Eigen::VectorXf vecf = boost::any_cast<Eigen::VectorXf const&>(alg->Zero(typeid(Eigen::VectorXf), 13));
  const Eigen::VectorXi veci = boost::any_cast<Eigen::VectorXi const&>(alg->Zero(typeid(Eigen::VectorXi), 13));
  EXPECT_DOUBLE_EQ(vecd.norm(), 0.0);
  EXPECT_EQ(vecd.size(), 13);
  EXPECT_DOUBLE_EQ(vecf.norm(), 0.0);
  EXPECT_EQ(vecf.size(), 13);
  EXPECT_DOUBLE_EQ(veci.norm(), 0.0);
  EXPECT_EQ(veci.size(), 13);

  // test the Eigen::Matrix zero
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix2d const&>(alg->Zero(typeid(Eigen::Matrix2d))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix2f const&>(alg->Zero(typeid(Eigen::Matrix2f))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix2i const&>(alg->Zero(typeid(Eigen::Matrix2i))).norm(), 0.0);

  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix3d const&>(alg->Zero(typeid(Eigen::Matrix3d))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix3f const&>(alg->Zero(typeid(Eigen::Matrix3f))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix3i const&>(alg->Zero(typeid(Eigen::Matrix3i))).norm(), 0.0);

  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix4d const&>(alg->Zero(typeid(Eigen::Matrix4d))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix4f const&>(alg->Zero(typeid(Eigen::Matrix4f))).norm(), 0.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<Eigen::Matrix4i const&>(alg->Zero(typeid(Eigen::Matrix4i))).norm(), 0.0);

  const Eigen::MatrixXd matd = boost::any_cast<Eigen::MatrixXd const&>(alg->Zero(typeid(Eigen::MatrixXd), 13, 5));
  const Eigen::MatrixXf matf = boost::any_cast<Eigen::MatrixXf const&>(alg->Zero(typeid(Eigen::MatrixXf), 13, 5));
  const Eigen::MatrixXi mati = boost::any_cast<Eigen::MatrixXi const&>(alg->Zero(typeid(Eigen::MatrixXi), 13, 5));
  EXPECT_DOUBLE_EQ(matd.norm(), 0.0);
  EXPECT_EQ(matd.rows(), 13);
  EXPECT_EQ(matd.cols(), 5);
  EXPECT_DOUBLE_EQ(matf.norm(), 0.0);
  EXPECT_EQ(matf.rows(), 13);
  EXPECT_EQ(matf.cols(), 5);
  EXPECT_DOUBLE_EQ(mati.norm(), 0.0);
  EXPECT_EQ(mati.rows(), 13);
  EXPECT_EQ(mati.cols(), 5);
}

TEST(AnyAlgebraTests, IsZero) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar zero
  EXPECT_TRUE(alg->IsZero((double)0.0));
  EXPECT_FALSE(alg->IsZero((double)1.0));
  EXPECT_TRUE(alg->IsZero((float)0.0));
  EXPECT_FALSE(alg->IsZero((float)-2.3));
  EXPECT_TRUE(alg->IsZero((int)0));
  EXPECT_FALSE(alg->IsZero((int)-6));
  EXPECT_TRUE(alg->IsZero((unsigned int)0));
  EXPECT_FALSE(alg->IsZero((unsigned int)18));

  EXPECT_TRUE(alg->IsZero((Eigen::Vector2d)Eigen::Vector2d::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector2d)Eigen::Vector2d::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Vector2f)Eigen::Vector2f::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector2f)Eigen::Vector2f::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Vector2i)Eigen::Vector2i::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector2i)Eigen::Vector2i::Random()));

  EXPECT_TRUE(alg->IsZero((Eigen::Vector3d)Eigen::Vector3d::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector3d)Eigen::Vector3d::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Vector3f)Eigen::Vector3f::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector3f)Eigen::Vector3f::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Vector3i)Eigen::Vector3i::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector3i)Eigen::Vector3i::Random()));

  EXPECT_TRUE(alg->IsZero((Eigen::Vector4d)Eigen::Vector4d::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector4d)Eigen::Vector4d::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Vector4f)Eigen::Vector4f::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector4f)Eigen::Vector4f::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Vector4i)Eigen::Vector4i::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Vector4i)Eigen::Vector4i::Random()));

  EXPECT_TRUE(alg->IsZero((Eigen::VectorXd)Eigen::VectorXd::Zero(21)));
  EXPECT_FALSE(alg->IsZero((Eigen::VectorXd)Eigen::VectorXd::Random(38)));
  EXPECT_TRUE(alg->IsZero((Eigen::VectorXf)Eigen::VectorXf::Zero(19)));
  EXPECT_FALSE(alg->IsZero((Eigen::VectorXf)Eigen::VectorXf::Random(22)));
  EXPECT_TRUE(alg->IsZero((Eigen::VectorXi)Eigen::VectorXi::Zero(36)));
  EXPECT_FALSE(alg->IsZero((Eigen::VectorXi)Eigen::VectorXi::Random(5)));

  EXPECT_TRUE(alg->IsZero((Eigen::Matrix2d)Eigen::Matrix2d::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix2d)Eigen::Matrix2d::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Matrix2f)Eigen::Matrix2f::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix2f)Eigen::Matrix2f::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Matrix2i)Eigen::Matrix2i::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix2i)Eigen::Matrix2i::Random()));

  EXPECT_TRUE(alg->IsZero((Eigen::Matrix3d)Eigen::Matrix3d::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix3d)Eigen::Matrix3d::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Matrix3f)Eigen::Matrix3f::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix3f)Eigen::Matrix3f::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Matrix3i)Eigen::Matrix3i::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix3i)Eigen::Matrix3i::Random()));

  EXPECT_TRUE(alg->IsZero((Eigen::Matrix4d)Eigen::Matrix4d::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix4d)Eigen::Matrix4d::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Matrix4f)Eigen::Matrix4f::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix4f)Eigen::Matrix4f::Random()));
  EXPECT_TRUE(alg->IsZero((Eigen::Matrix4i)Eigen::Matrix4i::Zero()));
  EXPECT_FALSE(alg->IsZero((Eigen::Matrix4i)Eigen::Matrix4i::Random()));

  EXPECT_TRUE(alg->IsZero((Eigen::MatrixXd)Eigen::MatrixXd::Zero(21, 5)));
  EXPECT_FALSE(alg->IsZero((Eigen::MatrixXd)Eigen::MatrixXd::Random(38, 13)));
  EXPECT_TRUE(alg->IsZero((Eigen::MatrixXf)Eigen::MatrixXf::Zero(19, 58)));
  EXPECT_FALSE(alg->IsZero((Eigen::MatrixXf)Eigen::MatrixXf::Random(22, 87)));
  EXPECT_TRUE(alg->IsZero((Eigen::MatrixXi)Eigen::MatrixXi::Zero(36, 23)));
  EXPECT_FALSE(alg->IsZero((Eigen::MatrixXi)Eigen::MatrixXi::Random(5, 9)));
}

TEST(AnyAlgebraTests, Identity) {
    auto alg = std::shared_ptr<AnyAlgebra>();

    // scalars
    EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Identity(typeid(double))), (double)1.0);
    EXPECT_DOUBLE_EQ(boost::any_cast<float const>(alg->Identity(typeid(float))), (float)1.0);
    EXPECT_DOUBLE_EQ(boost::any_cast<int const>(alg->Identity(typeid(int))), (int)1);
    EXPECT_DOUBLE_EQ(boost::any_cast<unsigned int const>(alg->Identity(typeid(unsigned int))), (unsigned int)1);
    
    // Eigen::Vector
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d>(alg->Identity(typeid(Eigen::Vector2d))).array()==Eigen::Matrix2d::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f>(alg->Identity(typeid(Eigen::Vector2f))).array()==Eigen::Matrix2f::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i>(alg->Identity(typeid(Eigen::Vector2i))).array()==Eigen::Matrix2i::Identity().array()).all());

    EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d>(alg->Identity(typeid(Eigen::Vector3d))).array()==Eigen::Matrix3d::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f>(alg->Identity(typeid(Eigen::Vector3f))).array()==Eigen::Matrix3f::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i>(alg->Identity(typeid(Eigen::Vector3i))).array()==Eigen::Matrix3i::Identity().array()).all());

    EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d>(alg->Identity(typeid(Eigen::Vector4d))).array()==Eigen::Matrix4d::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f>(alg->Identity(typeid(Eigen::Vector4f))).array()==Eigen::Matrix4f::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i>(alg->Identity(typeid(Eigen::Vector4i))).array()==Eigen::Matrix4i::Identity().array()).all());

    EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd>(alg->Identity(typeid(Eigen::VectorXd), 4, 9)).array()==Eigen::MatrixXd::Identity(4, 9).array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf>(alg->Identity(typeid(Eigen::VectorXf), 43, 21)).array()==Eigen::MatrixXf::Identity(43, 21).array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi>(alg->Identity(typeid(Eigen::VectorXi), 23, 23)).array()==Eigen::MatrixXi::Identity(23, 23).array()).all());

    // Eigen::Matrix
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d>(alg->Identity(typeid(Eigen::Matrix2d))).array()==Eigen::Matrix2d::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f>(alg->Identity(typeid(Eigen::Matrix2f))).array()==Eigen::Matrix2f::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i>(alg->Identity(typeid(Eigen::Matrix2i))).array()==Eigen::Matrix2i::Identity().array()).all());

    EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d>(alg->Identity(typeid(Eigen::Matrix3d))).array()==Eigen::Matrix3d::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f>(alg->Identity(typeid(Eigen::Matrix3f))).array()==Eigen::Matrix3f::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i>(alg->Identity(typeid(Eigen::Matrix3i))).array()==Eigen::Matrix3i::Identity().array()).all());

    EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d>(alg->Identity(typeid(Eigen::Matrix4d))).array()==Eigen::Matrix4d::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f>(alg->Identity(typeid(Eigen::Matrix4f))).array()==Eigen::Matrix4f::Identity().array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i>(alg->Identity(typeid(Eigen::Matrix4i))).array()==Eigen::Matrix4i::Identity().array()).all());

    EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd>(alg->Identity(typeid(Eigen::MatrixXd), 4, 9)).array()==Eigen::MatrixXd::Identity(4, 9).array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf>(alg->Identity(typeid(Eigen::MatrixXf), 43, 21)).array()==Eigen::MatrixXf::Identity(43, 21).array()).all());
    EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi>(alg->Identity(typeid(Eigen::MatrixXi), 23, 23)).array()==Eigen::MatrixXi::Identity(23, 23).array()).all());
}

TEST(AnyAlgebraTests, Norm) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test scalar norm
  EXPECT_DOUBLE_EQ(alg->Norm(-4.0), 4.0);

  // test Eigen::Vector norm
  const Eigen::Vector2d vec2(2.0, 2.0);
  EXPECT_DOUBLE_EQ(alg->Norm(vec2), vec2.norm());
  const Eigen::Vector3d vec3(2.0, 2.0, 2.0);
  EXPECT_DOUBLE_EQ(alg->Norm(vec3), vec3.norm());
  const Eigen::Vector4d vec4(2.0, 2.0, 2.0, 2.0);
  EXPECT_DOUBLE_EQ(alg->Norm(vec4), vec4.norm());
  const Eigen::VectorXd vec = Eigen::VectorXd::Random(8);
  EXPECT_DOUBLE_EQ(alg->Norm(vec), vec.norm());

  // test Eigen::Matrix norm
  const Eigen::Matrix2d mat2 = Eigen::Matrix2d::Random();
  EXPECT_DOUBLE_EQ(alg->Norm(mat2), mat2.norm());
  const Eigen::Matrix3d mat3 = Eigen::Matrix3d::Random();
  EXPECT_DOUBLE_EQ(alg->Norm(mat3), mat3.norm());
  const Eigen::Matrix4d mat4 = Eigen::Matrix4d::Random();
  EXPECT_DOUBLE_EQ(alg->Norm(mat4), mat4.norm());
  const Eigen::MatrixXd mat = Eigen::MatrixXd::Random(8, 4);
  EXPECT_DOUBLE_EQ(alg->Norm(mat), mat.norm());
}

TEST(AnyAlgebraTests, InnerProduct) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;
  EXPECT_EQ(alg->InnerProduct(xd, xd), xd*xd);
  EXPECT_EQ(alg->InnerProduct(xd, xf), xd*xf);
  EXPECT_EQ(alg->InnerProduct(xd, xi), xd*xi);
  EXPECT_EQ(alg->InnerProduct(xd, xui), xd*xui);
  EXPECT_EQ(alg->InnerProduct(xf, xd), xf*xd);
  EXPECT_EQ(alg->InnerProduct(xf, xf), xf*xf);
  EXPECT_EQ(alg->InnerProduct(xf, xi), xf*xi);
  EXPECT_EQ(alg->InnerProduct(xf, xui), xf*xui);
  EXPECT_EQ(alg->InnerProduct(xi, xd), xi*xd);
  EXPECT_EQ(alg->InnerProduct(xi, xf), xi*xf);
  EXPECT_EQ(alg->InnerProduct(xi, xi), xi*xi);
  EXPECT_EQ(alg->InnerProduct(xi, xui), xi*(double)xui);
  EXPECT_EQ(alg->InnerProduct(xui, xd), xui*xd);
  EXPECT_EQ(alg->InnerProduct(xui, xf), xui*xf);
  EXPECT_EQ(alg->InnerProduct(xui, xi), (double)xui*xi);
  EXPECT_EQ(alg->InnerProduct(xui, xui), xui*xui);

  const Eigen::Vector2d vec2(2.0, 3.0);
  const Eigen::Vector3d vec3(2.0, 3.0, 4.0);
  const Eigen::Vector4d vec4(2.0, 3.0, 4.0, 5.0);
  const Eigen::VectorXd vec = Eigen::VectorXd::Random(3);

  EXPECT_EQ(alg->InnerProduct(vec2, vec2), vec2.dot(vec2));
  EXPECT_EQ(alg->InnerProduct(vec3, vec3), vec3.dot(vec3));
  EXPECT_EQ(alg->InnerProduct(vec3, vec), vec3.dot(vec));
  EXPECT_EQ(alg->InnerProduct(vec4, vec4), vec4.dot(vec4));
  EXPECT_EQ(alg->InnerProduct(vec, vec), vec.dot(vec));
}

TEST(AnyAlgebraTests, Add) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Add(xd, xd)), xd+xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Add(xd, xf)), xd+xf);
  EXPECT_EQ(boost::any_cast<double const>(alg->Add(xd, xi)), xd+xi);
  EXPECT_EQ(boost::any_cast<double const>(alg->Add(xd, xui)), xd+xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Add(xf, xd)), xf+xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Add(xf, xf)), xf+xf);
  EXPECT_EQ(boost::any_cast<float const>(alg->Add(xf, xi)), xf+xi);
  EXPECT_EQ(boost::any_cast<float const>(alg->Add(xf, xui)), xf+xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Add(xi, xd)), xi+xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Add(xi, xf)), xi+xf);
  EXPECT_EQ(boost::any_cast<int const>(alg->Add(xi, xi)), xi+xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Add(xi, xui)), xi+xui); 
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Add(xui, xd)), xui+xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Add(xui, xf)), xui+xf);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Add(xui, xi)), xui+xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Add(xui, xui)), xui+xui);

  // test the Eigen::Vectors
  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random();
  const Eigen::VectorXd vec2Xd = Eigen::VectorXd::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d const&>(alg->Add(vec2d, vec2d)).array()==(vec2d+vec2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d const&>(alg->Add(vec2d, vec2Xd)).array()==(vec2d+vec2Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Add(vec2Xd, vec2d)).array()==(vec2d+vec2Xd).array()).all());

  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random();
  const Eigen::VectorXf vec2Xf = Eigen::VectorXf::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f const&>(alg->Add(vec2f, vec2f)).array()==(vec2f+vec2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f const&>(alg->Add(vec2f, vec2Xf)).array()==(vec2f+vec2Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Add(vec2Xf, vec2f)).array()==(vec2f+vec2Xf).array()).all());

  const Eigen::Vector2i vec2i = Eigen::Vector2i::Random();
  const Eigen::VectorXi vec2Xi = Eigen::VectorXi::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Add(vec2i, vec2i)).array()==(vec2i+vec2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Add(vec2i, vec2Xi)).array()==(vec2i+vec2Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Add(vec2Xi, vec2i)).array()==(vec2i+vec2Xi).array()).all());
  
  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random();
  const Eigen::VectorXd vec3Xd = Eigen::VectorXd::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d const&>(alg->Add(vec3d, vec3d)).array()==(vec3d+vec3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d const&>(alg->Add(vec3d, vec3Xd)).array()==(vec3d+vec3Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Add(vec3Xd, vec3d)).array()==(vec3d+vec3Xd).array()).all());

  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random();
  const Eigen::VectorXf vec3Xf = Eigen::VectorXf::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f const&>(alg->Add(vec3f, vec3f)).array()==(vec3f+vec3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f const&>(alg->Add(vec3f, vec3Xf)).array()==(vec3f+vec3Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Add(vec3Xf, vec3f)).array()==(vec3f+vec3Xf).array()).all());

  const Eigen::Vector3i vec3i = Eigen::Vector3i::Random();
  const Eigen::VectorXi vec3Xi = Eigen::VectorXi::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Add(vec3i, vec3i)).array()==(vec3i+vec3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Add(vec3i, vec3Xi)).array()==(vec3i+vec3Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Add(vec3Xi, vec3i)).array()==(vec3i+vec3Xi).array()).all());

  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random();
  const Eigen::VectorXd vec4Xd = Eigen::VectorXd::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d const&>(alg->Add(vec4d, vec4d)).array()==(vec4d+vec4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d const&>(alg->Add(vec4d, vec4Xd)).array()==(vec4d+vec4Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Add(vec4Xd, vec4d)).array()==(vec4d+vec4Xd).array()).all());

  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random();
  const Eigen::VectorXf vec4Xf = Eigen::VectorXf::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f const&>(alg->Add(vec4f, vec4f)).array()==(vec4f+vec4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f const&>(alg->Add(vec4f, vec4Xf)).array()==(vec4f+vec4Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Add(vec4Xf, vec4f)).array()==(vec4f+vec4Xf).array()).all());

  const Eigen::Vector4i vec4i = Eigen::Vector4i::Random();
  const Eigen::VectorXi vec4Xi = Eigen::VectorXi::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Add(vec4i, vec4i)).array()==(vec4i+vec4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Add(vec4i, vec4Xi)).array()==(vec4i+vec4Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Add(vec4Xi, vec4i)).array()==(vec4i+vec4Xi).array()).all());

  const Eigen::VectorXd vecd = Eigen::VectorXd::Random(8);
  const Eigen::VectorXf vecf = Eigen::VectorXf::Random(9);
  const Eigen::VectorXi veci = Eigen::VectorXi::Random(5);
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Add(vecd, vecd)).array()==(vecd+vecd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Add(vecf, vecf)).array()==(vecf+vecf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Add(veci, veci)).array()==(veci+veci).array()).all());

  // test the Eigen::Matrices
  const Eigen::Matrix2d mat2d = Eigen::Matrix2d::Random();
  const Eigen::MatrixXd mat2Xd = Eigen::MatrixXd::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Add(mat2d, mat2d)).array()==(mat2d+mat2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Add(mat2d, mat2Xd)).array()==(mat2d+mat2Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Add(mat2Xd, mat2d)).array()==(mat2d+mat2Xd).array()).all());

  const Eigen::Matrix2f mat2f = Eigen::Matrix2f::Random();
  const Eigen::MatrixXf mat2Xf = Eigen::MatrixXf::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Add(mat2f, mat2f)).array()==(mat2f+mat2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Add(mat2f, mat2Xf)).array()==(mat2f+mat2Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Add(mat2Xf, mat2f)).array()==(mat2f+mat2Xf).array()).all());

  const Eigen::Matrix2i mat2i = Eigen::Matrix2i::Random();
  const Eigen::MatrixXi mat2Xi = Eigen::MatrixXi::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Add(mat2i, mat2i)).array()==(mat2i+mat2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Add(mat2i, mat2Xi)).array()==(mat2i+mat2Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Add(mat2Xi, mat2i)).array()==(mat2i+mat2Xi).array()).all());
  
  const Eigen::Matrix3d mat3d = Eigen::Matrix3d::Random();
  const Eigen::MatrixXd mat3Xd = Eigen::MatrixXd::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Add(mat3d, mat3d)).array()==(mat3d+mat3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Add(mat3d, mat3Xd)).array()==(mat3d+mat3Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Add(mat3Xd, mat3d)).array()==(mat3d+mat3Xd).array()).all());

  const Eigen::Matrix3f mat3f = Eigen::Matrix3f::Random();
  const Eigen::MatrixXf mat3Xf = Eigen::MatrixXf::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Add(mat3f, mat3f)).array()==(mat3f+mat3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Add(mat3f, mat3Xf)).array()==(mat3f+mat3Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Add(mat3Xf, mat3f)).array()==(mat3f+mat3Xf).array()).all());

  const Eigen::Matrix3i mat3i = Eigen::Matrix3i::Random();
  const Eigen::MatrixXi mat3Xi = Eigen::MatrixXi::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Add(mat3i, mat3i)).array()==(mat3i+mat3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Add(mat3i, mat3Xi)).array()==(mat3i+mat3Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Add(mat3Xi, mat3i)).array()==(mat3i+mat3Xi).array()).all());

  const Eigen::Matrix4d mat4d = Eigen::Matrix4d::Random();
  const Eigen::MatrixXd mat4Xd = Eigen::MatrixXd::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Add(mat4d, mat4d)).array()==(mat4d+mat4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Add(mat4d, mat4Xd)).array()==(mat4d+mat4Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Add(mat4Xd, mat4d)).array()==(mat4d+mat4Xd).array()).all());

  const Eigen::Matrix4f mat4f = Eigen::Matrix4f::Random();
  const Eigen::MatrixXf mat4Xf = Eigen::MatrixXf::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Add(mat4f, mat4f)).array()==(mat4f+mat4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Add(mat4f, mat4Xf)).array()==(mat4f+mat4Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Add(mat4Xf, mat4f)).array()==(mat4f+mat4Xf).array()).all());

  const Eigen::Matrix4i mat4i = Eigen::Matrix4i::Random();
  const Eigen::MatrixXi mat4Xi = Eigen::MatrixXi::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Add(mat4i, mat4i)).array()==(mat4i+mat4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Add(mat4i, mat4Xi)).array()==(mat4i+mat4Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Add(mat4Xi, mat4i)).array()==(mat4i+mat4Xi).array()).all());

  const Eigen::MatrixXd matd = Eigen::MatrixXd::Random(8,62);
  const Eigen::MatrixXf matf = Eigen::MatrixXf::Random(9,2);
  const Eigen::MatrixXi mati = Eigen::MatrixXi::Random(5,5);
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Add(matd, matd)).array()==(matd+matd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Add(matf, matf)).array()==(matf+matf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Add(mati, mati)).array()==(mati+mati).array()).all());
}

TEST(AnyAlgebraTests, Subtract) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Subtract(xd, xd)), xd-xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Subtract(xd, xf)), xd-xf);
  EXPECT_EQ(boost::any_cast<double const>(alg->Subtract(xd, xi)), xd-xi);
  EXPECT_EQ(boost::any_cast<double const>(alg->Subtract(xd, xui)), xd-xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Subtract(xf, xd)), xf-xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Subtract(xf, xf)), xf-xf);
  EXPECT_EQ(boost::any_cast<float const>(alg->Subtract(xf, xi)), xf-xi);
  EXPECT_EQ(boost::any_cast<float const>(alg->Subtract(xf, xui)), xf-xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Subtract(xi, xd)), xi-xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Subtract(xi, xf)), xi-xf);
  EXPECT_EQ(boost::any_cast<int const>(alg->Subtract(xi, xi)), xi-xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Subtract(xi, xui)), xi-xui); 
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Subtract(xui, xd)), xui-xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Subtract(xui, xf)), xui-xf);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Subtract(xui, xi)), xui-xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Subtract(xui, xui)), xui-xui);

  // test the Eigen::Vectors
  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random();
  const Eigen::VectorXd vec2Xd = Eigen::VectorXd::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d const&>(alg->Subtract(vec2d, vec2d)).array()==(vec2d-vec2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d const&>(alg->Subtract(vec2d, vec2Xd)).array()==(vec2d-vec2Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Subtract(vec2Xd, vec2d)).array()==(vec2Xd-vec2d).array()).all());

  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random();
  const Eigen::VectorXf vec2Xf = Eigen::VectorXf::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f const&>(alg->Subtract(vec2f, vec2f)).array()==(vec2f-vec2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f const&>(alg->Subtract(vec2f, vec2Xf)).array()==(vec2f-vec2Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Subtract(vec2Xf, vec2f)).array()==(vec2Xf-vec2f).array()).all());

  const Eigen::Vector2i vec2i = Eigen::Vector2i::Random();
  const Eigen::VectorXi vec2Xi = Eigen::VectorXi::Random(2);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Subtract(vec2i, vec2i)).array()==(vec2i-vec2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Subtract(vec2i, vec2Xi)).array()==(vec2i-vec2Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Subtract(vec2Xi, vec2i)).array()==(vec2Xi-vec2i).array()).all());
  
  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random();
  const Eigen::VectorXd vec3Xd = Eigen::VectorXd::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d const&>(alg->Subtract(vec3d, vec3d)).array()==(vec3d-vec3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d const&>(alg->Subtract(vec3d, vec3Xd)).array()==(vec3d-vec3Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Subtract(vec3Xd, vec3d)).array()==(vec3Xd-vec3d).array()).all());

  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random();
  const Eigen::VectorXf vec3Xf = Eigen::VectorXf::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f const&>(alg->Subtract(vec3f, vec3f)).array()==(vec3f-vec3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f const&>(alg->Subtract(vec3f, vec3Xf)).array()==(vec3f-vec3Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Subtract(vec3Xf, vec3f)).array()==(vec3Xf-vec3f).array()).all());

  const Eigen::Vector3i vec3i = Eigen::Vector3i::Random();
  const Eigen::VectorXi vec3Xi = Eigen::VectorXi::Random(3);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Subtract(vec3i, vec3i)).array()==(vec3i-vec3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Subtract(vec3i, vec3Xi)).array()==(vec3i-vec3Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Subtract(vec3Xi, vec3i)).array()==(vec3Xi-vec3i).array()).all());

  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random();
  const Eigen::VectorXd vec4Xd = Eigen::VectorXd::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d const&>(alg->Subtract(vec4d, vec4d)).array()==(vec4d-vec4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d const&>(alg->Subtract(vec4d, vec4Xd)).array()==(vec4d-vec4Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Subtract(vec4Xd, vec4d)).array()==(vec4Xd-vec4d).array()).all());

  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random();
  const Eigen::VectorXf vec4Xf = Eigen::VectorXf::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f const&>(alg->Subtract(vec4f, vec4f)).array()==(vec4f-vec4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f const&>(alg->Subtract(vec4f, vec4Xf)).array()==(vec4f-vec4Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Subtract(vec4Xf, vec4f)).array()==(vec4Xf-vec4f).array()).all());

  const Eigen::Vector4i vec4i = Eigen::Vector4i::Random();
  const Eigen::VectorXi vec4Xi = Eigen::VectorXi::Random(4);
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Subtract(vec4i, vec4i)).array()==(vec4i-vec4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Subtract(vec4i, vec4Xi)).array()==(vec4i-vec4Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Subtract(vec4Xi, vec4i)).array()==(vec4Xi-vec4i).array()).all());

  const Eigen::VectorXd vecd = Eigen::VectorXd::Random(8);
  const Eigen::VectorXf vecf = Eigen::VectorXf::Random(9);
  const Eigen::VectorXi veci = Eigen::VectorXi::Random(5);
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Subtract(vecd, vecd)).array()==(vecd-vecd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Subtract(vecf, vecf)).array()==(vecf-vecf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Subtract(veci, veci)).array()==(veci-veci).array()).all());

  // test the Eigen::Matrices
  const Eigen::Matrix2d mat2d = Eigen::Matrix2d::Random();
  const Eigen::MatrixXd mat2Xd = Eigen::MatrixXd::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Subtract(mat2d, mat2d)).array()==(mat2d-mat2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Subtract(mat2d, mat2Xd)).array()==(mat2d-mat2Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Subtract(mat2Xd, mat2d)).array()==(mat2Xd-mat2d).array()).all());

  const Eigen::Matrix2f mat2f = Eigen::Matrix2f::Random();
  const Eigen::MatrixXf mat2Xf = Eigen::MatrixXf::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Subtract(mat2f, mat2f)).array()==(mat2f-mat2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Subtract(mat2f, mat2Xf)).array()==(mat2f-mat2Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Subtract(mat2Xf, mat2f)).array()==(mat2Xf-mat2f).array()).all());

  const Eigen::Matrix2i mat2i = Eigen::Matrix2i::Random();
  const Eigen::MatrixXi mat2Xi = Eigen::MatrixXi::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Subtract(mat2i, mat2i)).array()==(mat2i-mat2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Subtract(mat2i, mat2Xi)).array()==(mat2i-mat2Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Subtract(mat2Xi, mat2i)).array()==(mat2Xi-mat2i).array()).all());
  
  const Eigen::Matrix3d mat3d = Eigen::Matrix3d::Random();
  const Eigen::MatrixXd mat3Xd = Eigen::MatrixXd::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Subtract(mat3d, mat3d)).array()==(mat3d-mat3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Subtract(mat3d, mat3Xd)).array()==(mat3d-mat3Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Subtract(mat3Xd, mat3d)).array()==(mat3Xd-mat3d).array()).all());

  const Eigen::Matrix3f mat3f = Eigen::Matrix3f::Random();
  const Eigen::MatrixXf mat3Xf = Eigen::MatrixXf::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Subtract(mat3f, mat3f)).array()==(mat3f-mat3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Subtract(mat3f, mat3Xf)).array()==(mat3f-mat3Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Subtract(mat3Xf, mat3f)).array()==(mat3Xf-mat3f).array()).all());

  const Eigen::Matrix3i mat3i = Eigen::Matrix3i::Random();
  const Eigen::MatrixXi mat3Xi = Eigen::MatrixXi::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Subtract(mat3i, mat3i)).array()==(mat3i-mat3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Subtract(mat3i, mat3Xi)).array()==(mat3i-mat3Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Subtract(mat3Xi, mat3i)).array()==(mat3Xi-mat3i).array()).all());

  const Eigen::Matrix4d mat4d = Eigen::Matrix4d::Random();
  const Eigen::MatrixXd mat4Xd = Eigen::MatrixXd::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Subtract(mat4d, mat4d)).array()==(mat4d-mat4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Subtract(mat4d, mat4Xd)).array()==(mat4d-mat4Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Subtract(mat4Xd, mat4d)).array()==(mat4Xd-mat4d).array()).all());

  const Eigen::Matrix4f mat4f = Eigen::Matrix4f::Random();
  const Eigen::MatrixXf mat4Xf = Eigen::MatrixXf::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Subtract(mat4f, mat4f)).array()==(mat4f-mat4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Subtract(mat4f, mat4Xf)).array()==(mat4f-mat4Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Subtract(mat4Xf, mat4f)).array()==(mat4Xf-mat4f).array()).all());

  const Eigen::Matrix4i mat4i = Eigen::Matrix4i::Random();
  const Eigen::MatrixXi mat4Xi = Eigen::MatrixXi::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Subtract(mat4i, mat4i)).array()==(mat4i-mat4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Subtract(mat4i, mat4Xi)).array()==(mat4i-mat4Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Subtract(mat4Xi, mat4i)).array()==(mat4Xi-mat4i).array()).all());

  const Eigen::MatrixXd matd = Eigen::MatrixXd::Random(8,62);
  const Eigen::MatrixXf matf = Eigen::MatrixXf::Random(9,2);
  const Eigen::MatrixXi mati = Eigen::MatrixXi::Random(5,5);
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Subtract(matd, matd)).array()==(matd-matd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Subtract(matf, matf)).array()==(matf-matf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Subtract(mati, mati)).array()==(mati-mati).array()).all());
}

TEST(AnyAlgebraTests, Multiply) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Multiply(xd, xd)), xd*xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Multiply(xd, xf)), xd*xf);
  EXPECT_EQ(boost::any_cast<double const>(alg->Multiply(xd, xi)), xd*xi);
  EXPECT_EQ(boost::any_cast<double const>(alg->Multiply(xd, xui)), xd*xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Multiply(xf, xd)), xf*xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Multiply(xf, xf)), xf*xf);
  EXPECT_EQ(boost::any_cast<float const>(alg->Multiply(xf, xi)), xf*xi);
  EXPECT_EQ(boost::any_cast<float const>(alg->Multiply(xf, xui)), xf*xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Multiply(xi, xd)), xi*xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Multiply(xi, xf)), xi*xf);
  EXPECT_EQ(boost::any_cast<int const>(alg->Multiply(xi, xi)), xi*xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Multiply(xi, xui)), xi*xui); 
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Multiply(xui, xd)), xui*xd);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Multiply(xui, xf)), xui*xf);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Multiply(xui, xi)), xui*xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Multiply(xui, xui)), xui*xui);

  // test the Eigen::Vectors
  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d const&>(alg->Multiply(xd, vec2d)).array()==(xd*vec2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d const&>(alg->Multiply(vec2d, xd)).array()==(xd*vec2d).array()).all());

  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f const&>(alg->Multiply(xf, vec2f)).array()==(xf*vec2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f const&>(alg->Multiply(vec2f, xf)).array()==(xf*vec2f).array()).all());

  const Eigen::Vector2i vec2i = Eigen::Vector2i::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Multiply(xi, vec2i)).array()==(xi*vec2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Multiply(vec2i, xi)).array()==(xi*vec2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Multiply(xui, vec2i)).array()==(xui*vec2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i const&>(alg->Multiply(vec2i, xui)).array()==(xui*vec2i).array()).all());
  
  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d const&>(alg->Multiply(xd, vec3d)).array()==(xd*vec3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d const&>(alg->Multiply(vec3d, xd)).array()==(xd*vec3d).array()).all());

  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f const&>(alg->Multiply(xf, vec3f)).array()==(xf*vec3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f const&>(alg->Multiply(vec3f, xf)).array()==(xf*vec3f).array()).all());

  const Eigen::Vector3i vec3i = Eigen::Vector3i::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Multiply(xi, vec3i)).array()==(xi*vec3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Multiply(vec3i, xi)).array()==(xi*vec3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Multiply(xui, vec3i)).array()==(xui*vec3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i const&>(alg->Multiply(vec3i, xui)).array()==(xui*vec3i).array()).all());

  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d const&>(alg->Multiply(xd, vec4d)).array()==(xd*vec4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d const&>(alg->Multiply(vec4d, xd)).array()==(xd*vec4d).array()).all());

  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f const&>(alg->Multiply(xf, vec4f)).array()==(xf*vec4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f const&>(alg->Multiply(vec4f, xf)).array()==(xf*vec4f).array()).all());

  const Eigen::Vector4i vec4i = Eigen::Vector4i::Random();
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Multiply(xi, vec4i)).array()==(xi*vec4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Multiply(vec4i, xi)).array()==(xi*vec4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Multiply(xui, vec4i)).array()==(xui*vec4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i const&>(alg->Multiply(vec4i, xui)).array()==(xui*vec4i).array()).all());

  const Eigen::VectorXd vecd = Eigen::VectorXd::Random(12);
  const Eigen::VectorXf vecf = Eigen::VectorXf::Random(9);
  const Eigen::VectorXi veci = Eigen::VectorXi::Random(5);
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Multiply(xd, vecd)).array()==(xd*vecd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd const&>(alg->Multiply(vecd, xd)).array()==(xd*vecd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Multiply(xf, vecf)).array()==(xf*vecf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf const&>(alg->Multiply(vecf, xf)).array()==(xf*vecf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Multiply(xi, veci)).array()==(xi*veci).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Multiply(veci, xi)).array()==(xi*veci).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Multiply(xui, veci)).array()==(xui*veci).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi const&>(alg->Multiply(veci, xui)).array()==(xui*veci).array()).all());

  // test the Eigen::Matrices
  const Eigen::Matrix2d mat2d = Eigen::Matrix2d::Random();
  const Eigen::MatrixXd mat2Xd = Eigen::MatrixXd::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Multiply(xd, mat2d)).array()==(xd*mat2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Multiply(mat2d, xd)).array()==(xd*mat2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Multiply(mat2d, mat2d)).array()==(mat2d*mat2d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2d const&>(alg->Multiply(mat2d, mat2Xd)).array()==(mat2d*mat2Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Multiply(mat2Xd, mat2d)).array()==(mat2Xd*mat2d).array()).all());

  const Eigen::Matrix2f mat2f = Eigen::Matrix2f::Random();
  const Eigen::MatrixXf mat2Xf = Eigen::MatrixXf::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Multiply(xf, mat2f)).array()==(xf*mat2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Multiply(mat2f, xf)).array()==(xf*mat2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Multiply(mat2f, mat2f)).array()==(mat2f*mat2f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2f const&>(alg->Multiply(mat2f, mat2Xf)).array()==(mat2f*mat2Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Multiply(mat2Xf, mat2f)).array()==(mat2Xf*mat2f).array()).all());

  const Eigen::Matrix2i mat2i = Eigen::Matrix2i::Random();
  const Eigen::MatrixXi mat2Xi = Eigen::MatrixXi::Random(2,2);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Multiply(xi, mat2i)).array()==(xi*mat2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Multiply(mat2i, xi)).array()==(xi*mat2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Multiply(xui, mat2i)).array()==(xui*mat2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Multiply(mat2i, xui)).array()==(xui*mat2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Multiply(mat2i, mat2i)).array()==(mat2i*mat2i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix2i const&>(alg->Multiply(mat2i, mat2Xi)).array()==(mat2i*mat2Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(mat2Xi, mat2i)).array()==(mat2Xi*mat2i).array()).all());
  
  const Eigen::Matrix3d mat3d = Eigen::Matrix3d::Random();
  const Eigen::MatrixXd mat3Xd = Eigen::MatrixXd::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Multiply(xd, mat3d)).array()==(xd*mat3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Multiply(mat3d, xd)).array()==(xd*mat3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Multiply(mat3d, mat3d)).array()==(mat3d*mat3d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3d const&>(alg->Multiply(mat3d, mat3Xd)).array()==(mat3d*mat3Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Multiply(mat3Xd, mat3d)).array()==(mat3Xd*mat3d).array()).all());

  const Eigen::Matrix3f mat3f = Eigen::Matrix3f::Random();
  const Eigen::MatrixXf mat3Xf = Eigen::MatrixXf::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Multiply(xf, mat3f)).array()==(xf*mat3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Multiply(mat3f, xf)).array()==(xf*mat3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Multiply(mat3f, mat3f)).array()==(mat3f*mat3f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3f const&>(alg->Multiply(mat3f, mat3Xf)).array()==(mat3f*mat3Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Multiply(mat3Xf, mat3f)).array()==(mat3Xf*mat3f).array()).all());

  const Eigen::Matrix3i mat3i = Eigen::Matrix3i::Random();
  const Eigen::MatrixXi mat3Xi = Eigen::MatrixXi::Random(3,3);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Multiply(xi, mat3i)).array()==(xi*mat3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Multiply(mat3i, xi)).array()==(xi*mat3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Multiply(xui, mat3i)).array()==(xui*mat3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Multiply(mat3i, xui)).array()==(xui*mat3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Multiply(mat3i, mat3i)).array()==(mat3i*mat3i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix3i const&>(alg->Multiply(mat3i, mat3Xi)).array()==(mat3i*mat3Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(mat3Xi, mat3i)).array()==(mat3Xi*mat3i).array()).all());

  const Eigen::Matrix4d mat4d = Eigen::Matrix4d::Random();
  const Eigen::MatrixXd mat4Xd = Eigen::MatrixXd::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Multiply(xd, mat4d)).array()==(xd*mat4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Multiply(mat4d, xd)).array()==(xd*mat4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Multiply(mat4d, mat4d)).array()==(mat4d*mat4d).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4d const&>(alg->Multiply(mat4d, mat4Xd)).array()==(mat4d*mat4Xd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Multiply(mat4Xd, mat4d)).array()==(mat4Xd*mat4d).array()).all());

  const Eigen::Matrix4f mat4f = Eigen::Matrix4f::Random();
  const Eigen::MatrixXf mat4Xf = Eigen::MatrixXf::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Multiply(xf, mat4f)).array()==(xf*mat4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Multiply(mat4f, xf)).array()==(xf*mat4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Multiply(mat4f, mat4f)).array()==(mat4f*mat4f).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4f const&>(alg->Multiply(mat4f, mat4Xf)).array()==(mat4f*mat4Xf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Multiply(mat4Xf, mat4f)).array()==(mat4Xf*mat4f).array()).all());

  const Eigen::Matrix4i mat4i = Eigen::Matrix4i::Random();
  const Eigen::MatrixXi mat4Xi = Eigen::MatrixXi::Random(4,4);
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Multiply(xi, mat4i)).array()==(xi*mat4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Multiply(mat4i, xi)).array()==(xi*mat4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Multiply(xui, mat4i)).array()==(xui*mat4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Multiply(mat4i, xui)).array()==(xui*mat4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Multiply(mat4i, mat4i)).array()==(mat4i*mat4i).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Matrix4i const&>(alg->Multiply(mat4i, mat4Xi)).array()==(mat4i*mat4Xi).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(mat4Xi, mat4i)).array()==(mat4Xi*mat4i).array()).all());

  const Eigen::MatrixXd matd = Eigen::MatrixXd::Random(8,8);
  const Eigen::MatrixXf matf = Eigen::MatrixXf::Random(9,9);
  const Eigen::MatrixXi mati = Eigen::MatrixXi::Random(5,5);
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Multiply(xd, matd)).array()==(xd*matd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Multiply(matd, xd)).array()==(xd*matd).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXd const&>(alg->Multiply(matd, matd)).array()==(matd*matd).array()).all());

  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Multiply(xf, matf)).array()==(xf*matf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Multiply(matf, xf)).array()==(xf*matf).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXf const&>(alg->Multiply(matf, matf)).array()==(matf*matf).array()).all());
  
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(xi, mati)).array()==(xi*mati).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(mati, xi)).array()==(xi*mati).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(xui, mati)).array()==(xui*mati).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(mati, xui)).array()==(xui*mati).array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::MatrixXi const&>(alg->Multiply(mati, mati)).array()==(mati*mati).array()).all());
}

TEST(AnyAlgebraTests, Apply) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;

  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Apply(xd, xd)), xd*xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Apply(xd, xf)), xf*xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Apply(xd, xi)), xi*xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Apply(xd, xui)), xui*xd);

  EXPECT_FLOAT_EQ(boost::any_cast<double const>(alg->Apply(xf, xd)), xd*xf); // use FLOAT_EQ because the precision is lower
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Apply(xf, xf)), xf*xf);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Apply(xf, xi)), xi*xf);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Apply(xf, xui)), xui*xf);
  
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Apply(xi, xd)), xd*xi);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Apply(xi, xf)), xf*xi);
  EXPECT_EQ(boost::any_cast<int const>(alg->Apply(xi, xi)), xi*xi);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Apply(xi, xui)), xui*xi);

  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->Apply(xui, xd)), xd*xui);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->Apply(xui, xf)), xf*xui);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Apply(xui, xi)), xi*xui);
  EXPECT_EQ(boost::any_cast<unsigned int const>(alg->Apply(xui, xui)), xui*xui);

  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random();
  const Eigen::VectorXd vec2Xd = Eigen::VectorXd::Random(2);
  const Eigen::Vector2d vm_2d2d = vec2d.asDiagonal()*vec2d;
  const Eigen::Vector2d vm_2d2Xd = vec2d.asDiagonal()*vec2Xd;
  const Eigen::Vector2d vm_2Xd2d = vec2Xd.asDiagonal()*vec2d;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d>(alg->Apply(vec2d, vec2d)).array()==vm_2d2d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->Apply(vec2d, vec2Xd)).array()==vm_2d2Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d>(alg->Apply(vec2Xd, vec2d)).array()==vm_2Xd2d.array()).all());

  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random();
  const Eigen::VectorXf vec2Xf = Eigen::VectorXf::Random(2);
  const Eigen::Vector2f vm_2f2f = vec2f.asDiagonal()*vec2f;
  const Eigen::Vector2f vm_2f2Xf = vec2f.asDiagonal()*vec2Xf;
  const Eigen::Vector2f vm_2Xf2f = vec2Xf.asDiagonal()*vec2f;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f>(alg->Apply(vec2f, vec2f)).array()==vm_2f2f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->Apply(vec2f, vec2Xf)).array()==vm_2f2Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f>(alg->Apply(vec2Xf, vec2f)).array()==vm_2Xf2f.array()).all());

  const Eigen::Vector2i vec2i = Eigen::Vector2i::Random();
  const Eigen::VectorXi vec2Xi = Eigen::VectorXi::Random(2);
  const Eigen::Vector2i vm_2i2i = vec2i.asDiagonal()*vec2i;
  const Eigen::Vector2i vm_2i2Xi = vec2i.asDiagonal()*vec2Xi;
  const Eigen::Vector2i vm_2Xi2i = vec2Xi.asDiagonal()*vec2i;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i>(alg->Apply(vec2i, vec2i)).array()==vm_2i2i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->Apply(vec2i, vec2Xi)).array()==vm_2i2Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i>(alg->Apply(vec2Xi, vec2i)).array()==vm_2Xi2i.array()).all());

  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random();
  const Eigen::VectorXd vec3Xd = Eigen::VectorXd::Random(3);
  const Eigen::Vector3d vm_3d3d = vec3d.asDiagonal()*vec3d;
  const Eigen::Vector3d vm_3d3Xd = vec3d.asDiagonal()*vec3Xd;
  const Eigen::Vector3d vm_3Xd3d = vec3Xd.asDiagonal()*vec3d;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d>(alg->Apply(vec3d, vec3d)).array()==vm_3d3d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->Apply(vec3d, vec3Xd)).array()==vm_3d3Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d>(alg->Apply(vec3Xd, vec3d)).array()==vm_3Xd3d.array()).all());

  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random();
  const Eigen::VectorXf vec3Xf = Eigen::VectorXf::Random(3);
  const Eigen::Vector3f vm_3f3f = vec3f.asDiagonal()*vec3f;
  const Eigen::Vector3f vm_3f3Xf = vec3f.asDiagonal()*vec3Xf;
  const Eigen::Vector3f vm_3Xf3f = vec3Xf.asDiagonal()*vec3f;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f>(alg->Apply(vec3f, vec3f)).array()==vm_3f3f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->Apply(vec3f, vec3Xf)).array()==vm_3f3Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f>(alg->Apply(vec3Xf, vec3f)).array()==vm_3Xf3f.array()).all());

  const Eigen::Vector3i vec3i = Eigen::Vector3i::Random();
  const Eigen::VectorXi vec3Xi = Eigen::VectorXi::Random(3);
  const Eigen::Vector3i vm_3i3i = vec3i.asDiagonal()*vec3i;
  const Eigen::Vector3i vm_3i3Xi = vec3i.asDiagonal()*vec3Xi;
  const Eigen::Vector3i vm_3Xi3i = vec3Xi.asDiagonal()*vec3i;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i>(alg->Apply(vec3i, vec3i)).array()==vm_3i3i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->Apply(vec3i, vec3Xi)).array()==vm_3i3Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i>(alg->Apply(vec3Xi, vec3i)).array()==vm_3Xi3i.array()).all());

  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random();
  const Eigen::VectorXd vec4Xd = Eigen::VectorXd::Random(4);
  const Eigen::Vector4d vm_4d4d = vec4d.asDiagonal()*vec4d;
  const Eigen::Vector4d vm_4d4Xd = vec4d.asDiagonal()*vec4Xd;
  const Eigen::Vector4d vm_4Xd4d = vec4Xd.asDiagonal()*vec4d;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d>(alg->Apply(vec4d, vec4d)).array()==vm_4d4d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->Apply(vec4d, vec4Xd)).array()==vm_4d4Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d>(alg->Apply(vec4Xd, vec4d)).array()==vm_4Xd4d.array()).all());

  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random();
  const Eigen::VectorXf vec4Xf = Eigen::VectorXf::Random(4);
  const Eigen::Vector4f vm_4f4f = vec4f.asDiagonal()*vec4f;
  const Eigen::Vector4f vm_4f4Xf = vec4f.asDiagonal()*vec4Xf;
  const Eigen::Vector4f vm_4Xf4f = vec4Xf.asDiagonal()*vec4f;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f>(alg->Apply(vec4f, vec4f)).array()==vm_4f4f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->Apply(vec4f, vec4Xf)).array()==vm_4f4Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f>(alg->Apply(vec4Xf, vec4f)).array()==vm_4Xf4f.array()).all());

  const Eigen::Vector4i vec4i = Eigen::Vector4i::Random();
  const Eigen::VectorXi vec4Xi = Eigen::VectorXi::Random(4);
  const Eigen::Vector4i vm_4i4i = vec4i.asDiagonal()*vec4i;
  const Eigen::Vector4i vm_4i4Xi = vec4i.asDiagonal()*vec4Xi;
  const Eigen::Vector4i vm_4Xi4i = vec4Xi.asDiagonal()*vec4i;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i>(alg->Apply(vec4i, vec4i)).array()==vm_4i4i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->Apply(vec4i, vec4Xi)).array()==vm_4i4Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i>(alg->Apply(vec4Xi, vec4i)).array()==vm_4Xi4i.array()).all());

  const Eigen::VectorXd vecd = Eigen::VectorXd::Random(13);
  const Eigen::VectorXd vm_dd = vecd.asDiagonal()*vecd;
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->Apply(vecd, vecd)).array()==vm_dd.array()).all());

  const Eigen::VectorXf vecf = Eigen::VectorXf::Random(13);
  const Eigen::VectorXf vm_ff = vecf.asDiagonal()*vecf;
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->Apply(vecf, vecf)).array()==vm_ff.array()).all());

  const Eigen::VectorXi veci = Eigen::VectorXi::Random(13);
  const Eigen::VectorXi vm_ii = veci.asDiagonal()*veci;
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->Apply(veci, veci)).array()==vm_ii.array()).all());
}

TEST(AnyAlgebraTests, ApplyInverse) {
  auto alg = std::shared_ptr<AnyAlgebra>();

  // test the scalar inner product
  double xd=4.0; float xf=-3.0; int xi=-2; unsigned int xui=8;

  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xd, xd)), 1.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xd, xf)), xf/xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xd, xi)), xi/xd);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xd, xui)), xui/xd);

  EXPECT_FLOAT_EQ(boost::any_cast<double const>(alg->ApplyInverse(xf, xd)), xd/xf); // use FLOAT_EQ because the precision is lower
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->ApplyInverse(xf, xf)), 1.0);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->ApplyInverse(xf, xi)), xi/xf);
  EXPECT_FLOAT_EQ(boost::any_cast<float const>(alg->ApplyInverse(xf, xui)), xui/xf);
  
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xi, xd)), xd/xi);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xi, xf)), xf/xi);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xi, xi)), 1.0);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xi, xui)), (double)xui/(double)xi);

  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xui, xd)), xd/xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xui, xf)), xf/xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xui, xi)), (double)xi/(double)xui);
  EXPECT_DOUBLE_EQ(boost::any_cast<double const>(alg->ApplyInverse(xui, xui)), 1);

  const Eigen::Vector2d vec2d = Eigen::Vector2d::Random();
  const Eigen::VectorXd vec2Xd = Eigen::VectorXd::Random(2);
  const Eigen::Vector2d vm_2d2d = (1.0/vec2d.array()).matrix().asDiagonal()*vec2d;
  const Eigen::Vector2d vm_2d2Xd = (1.0/vec2d.array()).matrix().asDiagonal()*vec2Xd;
  const Eigen::Vector2d vm_2Xd2d = (1.0/vec2Xd.array()).matrix().asDiagonal()*vec2d;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d>(alg->ApplyInverse(vec2d, vec2d)).array()==vm_2d2d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->ApplyInverse(vec2d, vec2Xd)).array()==vm_2d2Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2d>(alg->ApplyInverse(vec2Xd, vec2d)).array()==vm_2Xd2d.array()).all());

  const Eigen::Vector2f vec2f = Eigen::Vector2f::Random();
  const Eigen::VectorXf vec2Xf = Eigen::VectorXf::Random(2);
  const Eigen::Vector2f vm_2f2f = (1.0/vec2f.array()).matrix().asDiagonal()*vec2f;
  const Eigen::Vector2f vm_2f2Xf = (1.0/vec2f.array()).matrix().asDiagonal()*vec2Xf;
  const Eigen::Vector2f vm_2Xf2f = (1.0/vec2Xf.array()).matrix().asDiagonal()*vec2f;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f>(alg->ApplyInverse(vec2f, vec2f)).array()==vm_2f2f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->ApplyInverse(vec2f, vec2Xf)).array()==vm_2f2Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2f>(alg->ApplyInverse(vec2Xf, vec2f)).array()==vm_2Xf2f.array()).all());

  const Eigen::Vector2i vec2i = Eigen::Vector2i::Random();
  const Eigen::VectorXi vec2Xi = Eigen::VectorXi::Random(2);
  const Eigen::Vector2i vm_2i2i = (1/vec2i.array()).matrix().asDiagonal()*vec2i;
  const Eigen::Vector2i vm_2i2Xi = (1/vec2i.array()).matrix().asDiagonal()*vec2Xi;
  const Eigen::Vector2i vm_2Xi2i = (1/vec2Xi.array()).matrix().asDiagonal()*vec2i;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i>(alg->ApplyInverse(vec2i, vec2i)).array()==vm_2i2i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->ApplyInverse(vec2i, vec2Xi)).array()==vm_2i2Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector2i>(alg->ApplyInverse(vec2Xi, vec2i)).array()==vm_2Xi2i.array()).all());

  const Eigen::Vector3d vec3d = Eigen::Vector3d::Random();
  const Eigen::VectorXd vec3Xd = Eigen::VectorXd::Random(3);
  const Eigen::Vector3d vm_3d3d = (1.0/vec3d.array()).matrix().asDiagonal()*vec3d;
  const Eigen::Vector3d vm_3d3Xd = (1.0/vec3d.array()).matrix().asDiagonal()*vec3Xd;
  const Eigen::Vector3d vm_3Xd3d = (1.0/vec3Xd.array()).matrix().asDiagonal()*vec3d;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d>(alg->ApplyInverse(vec3d, vec3d)).array()==vm_3d3d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->ApplyInverse(vec3d, vec3Xd)).array()==vm_3d3Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3d>(alg->ApplyInverse(vec3Xd, vec3d)).array()==vm_3Xd3d.array()).all());

  const Eigen::Vector3f vec3f = Eigen::Vector3f::Random();
  const Eigen::VectorXf vec3Xf = Eigen::VectorXf::Random(3);
  const Eigen::Vector3f vm_3f3f = (1.0/vec3f.array()).matrix().asDiagonal()*vec3f;
  const Eigen::Vector3f vm_3f3Xf = (1.0/vec3f.array()).matrix().asDiagonal()*vec3Xf;
  const Eigen::Vector3f vm_3Xf3f = (1.0/vec3Xf.array()).matrix().asDiagonal()*vec3f;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f>(alg->ApplyInverse(vec3f, vec3f)).array()==vm_3f3f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->ApplyInverse(vec3f, vec3Xf)).array()==vm_3f3Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3f>(alg->ApplyInverse(vec3Xf, vec3f)).array()==vm_3Xf3f.array()).all());

  const Eigen::Vector3i vec3i = Eigen::Vector3i::Random();
  const Eigen::VectorXi vec3Xi = Eigen::VectorXi::Random(3);
  const Eigen::Vector3i vm_3i3i = (1/vec3i.array()).matrix().asDiagonal()*vec3i;
  const Eigen::Vector3i vm_3i3Xi = (1/vec3i.array()).matrix().asDiagonal()*vec3Xi;
  const Eigen::Vector3i vm_3Xi3i = (1/vec3Xi.array()).matrix().asDiagonal()*vec3i;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i>(alg->ApplyInverse(vec3i, vec3i)).array()==vm_3i3i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->ApplyInverse(vec3i, vec3Xi)).array()==vm_3i3Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector3i>(alg->ApplyInverse(vec3Xi, vec3i)).array()==vm_3Xi3i.array()).all());

  const Eigen::Vector4d vec4d = Eigen::Vector4d::Random();
  const Eigen::VectorXd vec4Xd = Eigen::VectorXd::Random(4);
  const Eigen::Vector4d vm_4d4d = (1.0/vec4d.array()).matrix().asDiagonal()*vec4d;
  const Eigen::Vector4d vm_4d4Xd = (1.0/vec4d.array()).matrix().asDiagonal()*vec4Xd;
  const Eigen::Vector4d vm_4Xd4d = (1.0/vec4Xd.array()).matrix().asDiagonal()*vec4d;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d>(alg->ApplyInverse(vec4d, vec4d)).array()==vm_4d4d.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->ApplyInverse(vec4d, vec4Xd)).array()==vm_4d4Xd.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4d>(alg->ApplyInverse(vec4Xd, vec4d)).array()==vm_4Xd4d.array()).all());

  const Eigen::Vector4f vec4f = Eigen::Vector4f::Random();
  const Eigen::VectorXf vec4Xf = Eigen::VectorXf::Random(4);
  const Eigen::Vector4f vm_4f4f = (1.0/vec4f.array()).matrix().asDiagonal()*vec4f;
  const Eigen::Vector4f vm_4f4Xf = (1.0/vec4f.array()).matrix().asDiagonal()*vec4Xf;
  const Eigen::Vector4f vm_4Xf4f = (1.0/vec4Xf.array()).matrix().asDiagonal()*vec4f;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f>(alg->ApplyInverse(vec4f, vec4f)).array()==vm_4f4f.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->ApplyInverse(vec4f, vec4Xf)).array()==vm_4f4Xf.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4f>(alg->ApplyInverse(vec4Xf, vec4f)).array()==vm_4Xf4f.array()).all());

  const Eigen::Vector4i vec4i = Eigen::Vector4i::Random();
  const Eigen::VectorXi vec4Xi = Eigen::VectorXi::Random(4);
  const Eigen::Vector4i vm_4i4i = (1/vec4i.array()).matrix().asDiagonal()*vec4i;
  const Eigen::Vector4i vm_4i4Xi = (1/vec4i.array()).matrix().asDiagonal()*vec4Xi;
  const Eigen::Vector4i vm_4Xi4i = (1/vec4Xi.array()).matrix().asDiagonal()*vec4i;
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i>(alg->ApplyInverse(vec4i, vec4i)).array()==vm_4i4i.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->ApplyInverse(vec4i, vec4Xi)).array()==vm_4i4Xi.array()).all());
  EXPECT_TRUE((boost::any_cast<Eigen::Vector4i>(alg->ApplyInverse(vec4Xi, vec4i)).array()==vm_4Xi4i.array()).all());

  const Eigen::VectorXd vecd = Eigen::VectorXd::Random(13);
  const Eigen::VectorXd vm_dd = (1.0/vecd.array()).matrix().asDiagonal()*vecd;
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXd>(alg->ApplyInverse(vecd, vecd)).array()==vm_dd.array()).all());

  const Eigen::VectorXf vecf = Eigen::VectorXf::Random(13);
  const Eigen::VectorXf vm_ff = (1.0/vecf.array()).matrix().asDiagonal()*vecf;
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXf>(alg->ApplyInverse(vecf, vecf)).array()==vm_ff.array()).all());

  const Eigen::VectorXi veci = Eigen::VectorXi::Random(13);
  const Eigen::VectorXi vm_ii = (1/veci.array()).matrix().asDiagonal()*veci;
  EXPECT_TRUE((boost::any_cast<Eigen::VectorXi>(alg->ApplyInverse(veci, veci)).array()==vm_ii.array()).all());
}
