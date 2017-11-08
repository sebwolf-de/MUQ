#include <gtest/gtest.h>

#include "MUQ/Utilities/LinearAlgebra/AnyAlgebra.h"

using namespace muq::Utilities;

TEST(AnyAlgebraTests, Size) {
  auto alg = std::shared_ptr<AnyAlgebra>();
    
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

TEST(AnyAlgebraTests, Add) {
  auto alg = std::shared_ptr<AnyAlgebra>();

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
