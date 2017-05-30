#include "gtest/gtest.h"

#include <Eigen/Core>

#include "MUQ/Modeling/RootfindingIVP.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;

/// The right hand side for a chemical kenetics ODE (from Sundials example)
class KineticsRHS : public WorkPiece {
public:

  /// Constructor
  inline KineticsRHS() : WorkPiece(std::vector<std::string>({typeid(N_Vector).name(), typeid(Eigen::Vector2d).name(), typeid(double).name()}), std::vector<std::string>(1, typeid(N_Vector).name())) {}

  inline virtual ~KineticsRHS() {}
  
private:

  inline virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override {
    outputs.resize(1);
    
    // dependent variables
    const N_Vector& y = boost::any_cast<const N_Vector>(inputs[0]);

    // parameters
    const Eigen::Vector2d& a = boost::any_cast<const Eigen::Vector2d>(inputs[1]);
    const double b = boost::any_cast<const double>(inputs[2]);

    // derivative wrt time
    outputs[0] = N_VNew_Serial(3);
    N_Vector& refout = boost::any_cast<N_Vector&>(outputs[0]);
    NV_Ith_S(refout, 0) = a(0)*NV_Ith_S(y, 0); // dy1dt
    NV_Ith_S(refout, 1) = a(1)*NV_Ith_S(y, 0)*NV_Ith_S(y, 1); // dy2dt
    NV_Ith_S(refout, 2) = b*b; // dy3dt
  }

  inline virtual void JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) override {
    // dependent variables
    const N_Vector& y = boost::any_cast<const N_Vector>(inputs[0]);

    // parameters
    const Eigen::Vector2d& a = boost::any_cast<const Eigen::Vector2d>(inputs[1]);
    const double b = boost::any_cast<const double>(inputs[2]);

    assert(wrtOut==0);
    if( wrtIn==0 ) { // dervative wrt y
      // get a reference to the jacobian matrix
      jacobian = NewDenseMat(3, 3);
      DlsMat& jac = boost::any_cast<DlsMat&>(*jacobian);
      SetToZero(jac);

      // implement the jacobian
      DENSE_ELEM(jac, 0, 0) = a(0);
      DENSE_ELEM(jac, 1, 0) = a(1)*NV_Ith_S(y, 1);
      DENSE_ELEM(jac, 1, 1) = a(1)*NV_Ith_S(y, 0);
    } else if( wrtIn==1 ) { // derivative wrt a
      // get a reference to the jacobian matrix
      jacobian = NewDenseMat(3, 2);
      DlsMat& jac = boost::any_cast<DlsMat&>(*jacobian);
      SetToZero(jac);

      DENSE_ELEM(jac, 0, 0) = NV_Ith_S(y, 0);
      DENSE_ELEM(jac, 1, 1) = NV_Ith_S(y, 0)*NV_Ith_S(y, 1);
    } else if( wrtIn==2 ) { // dervative wrt b
      // get a reference to the jacobian matrix
      jacobian = NewDenseMat(3, 1);
      DlsMat& jac = boost::any_cast<DlsMat&>(*jacobian);
      SetToZero(jac);

      DENSE_ELEM(jac, 2, 0) = 2.0*b;
    }
  }
};

/// We want to the integrate the kinetics probem until we find the root of this function (from a Sundials example)
class RootFunction : public WorkPiece {
public:

  /// Constructor
  inline RootFunction() : WorkPiece(std::vector<std::string>({typeid(N_Vector).name(), typeid(Eigen::Vector2d).name()}), std::vector<std::string>(2, typeid(double).name())) {}

  inline virtual ~RootFunction() {}
  
private:

  inline virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override {
    // dependent variables
    const N_Vector& y = boost::any_cast<const N_Vector>(inputs[0]);

    // parameters
    const Eigen::Vector2d& p = boost::any_cast<const Eigen::Vector2d>(inputs[1]);

    // compute the root function
    outputs.resize(2);
    outputs[0] = p(0); // root 1 (this will only be zero if p(0) is ...)
    outputs[1] = NV_Ith_S(y, 2)*NV_Ith_S(y, 2)-p(1); // root 2 
  }

  inline virtual void JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) override {
    /*// dependent variables
    const Eigen::Vector3d& y = boost::any_cast<const Eigen::Vector3d>(inputs[0]);

    // parameters
    const Eigen::Vector2d& p = boost::any_cast<const Eigen::Vector2d>(inputs[1]);

    // compute the jacobian 
    if( wrtIn==0 ) { // wrt the state y
      if( wrtOut==0 ) { // first output
	jacobian = (Eigen::MatrixXd)Eigen::MatrixXd::Zero(1, 3);
      } else if( wrtOut==1 ) { // second output
	jacobian = (Eigen::MatrixXd)Eigen::MatrixXd::Zero(1, 3);
	Eigen::MatrixXd& jac = boost::any_cast<Eigen::MatrixXd&>(*jacobian);

	jac(0, 2) = 2.0*y(2);
      }
    } else if( wrtIn==1 ) { // wrt the parameters p
      if( wrtOut==0 ) { // first output
	jacobian = (Eigen::MatrixXd)Eigen::MatrixXd::Zero(1, 2);
	Eigen::MatrixXd& jac = boost::any_cast<Eigen::MatrixXd&>(*jacobian);

	jac(0, 0) = 1.0;
      } else if( wrtOut==1 ) { // second output
	jacobian = (Eigen::MatrixXd)Eigen::MatrixXd::Zero(1, 2);
	Eigen::MatrixXd& jac = boost::any_cast<Eigen::MatrixXd&>(*jacobian);

	jac(0, 1) = -1.0;
      }
      }*/
  }
};

// test the RHS of the kinetics problem (not actually part of MUQ, but a good sanity check)
TEST(KineticsProblemTest, RHS) {
  // the right hand side of the ODE
  auto rhs = std::make_shared<KineticsRHS>();

  // inputs
  const Eigen::Vector3d y = Eigen::Vector3d::Random(3);
  N_Vector yvec = N_VNew_Serial(3);
  NV_Ith_S(yvec, 0) = y(0);
  NV_Ith_S(yvec, 1) = y(1);
  NV_Ith_S(yvec, 2) = y(2);
  const Eigen::Vector2d a = Eigen::Vector2d::Random(2);
  const double b = 3.5;

  { // test evaluate
    const std::vector<boost::any>& result = rhs->Evaluate(yvec, a, b);
    const N_Vector& rhs = boost::any_cast<const N_Vector>(result[0]);

    // check result
    EXPECT_DOUBLE_EQ(NV_Ith_S(rhs, 0), a(0)*y(0));
    EXPECT_DOUBLE_EQ(NV_Ith_S(rhs, 1), a(1)*y(0)*y(1));
    EXPECT_DOUBLE_EQ(NV_Ith_S(rhs, 2), b*b);
  }

  { // test jacobian
    // wrt to the state y
    const boost::any& jac0 = rhs->Jacobian(0, 0, yvec, a, b);
    const DlsMat& jac0ref = boost::any_cast<const DlsMat>(jac0);

    Eigen::MatrixXd expectedJac0 = Eigen::MatrixXd::Zero(3, 3);
    expectedJac0(0, 0) = a(0);
    expectedJac0(1, 0) = a(1)*y(1);
    expectedJac0(1, 1) = a(1)*y(0);

    EXPECT_EQ(jac0ref->M, 3); // check the number of rows
    EXPECT_EQ(jac0ref->N, 3); // check the number of cols

    // wrt to the parameters a
    const boost::any& jac1 = rhs->Jacobian(1, 0, yvec, a, b);
    const DlsMat& jac1ref = boost::any_cast<const DlsMat>(jac1);

    Eigen::MatrixXd expectedJac1 = Eigen::MatrixXd::Zero(3, 2);
    expectedJac1(0, 0) = y(0);
    expectedJac1(1, 1) = y(0)*y(1);

    EXPECT_EQ(jac1ref->M, 3); // check the number of rows
    EXPECT_EQ(jac1ref->N, 2); // check the number of cols

    // wrt to the parameters b
    const boost::any& jac2 = rhs->Jacobian(2, 0, yvec, a, b);
    const DlsMat& jac2ref = boost::any_cast<const DlsMat>(jac2);

    Eigen::MatrixXd expectedJac2 = Eigen::MatrixXd::Zero(3, 1);
    expectedJac2(2, 0) = 2.0*b;

    EXPECT_EQ(jac2ref->M, 3); // check the number of rows
    EXPECT_EQ(jac2ref->N, 1); // check the number of cols

    // check the values
    for( unsigned int i=0; i<3; ++i ) {
      for( unsigned int j=0; j<3; ++j ) {
	EXPECT_DOUBLE_EQ(DENSE_ELEM(jac0ref, i, j), expectedJac0(i,j));
      }

      for( unsigned int j=0; j<2; ++j ) {
	EXPECT_DOUBLE_EQ(DENSE_ELEM(jac1ref, i, j), expectedJac1(i,j));
      }

      EXPECT_DOUBLE_EQ(DENSE_ELEM(jac2ref, i, 0), expectedJac2(i,0));
    }
  }
}

// test the root function of the kinetics problem (not actually part of MUQ, but a good sanity check)
TEST(KineticsProblemTest, Root) {
  // we are trying to find the root of this function
  auto root = std::make_shared<RootFunction>();

  // inputs
  const Eigen::Vector3d y = Eigen::Vector3d::Random(3);
  N_Vector yvec = N_VNew_Serial(3);
  NV_Ith_S(yvec, 0) = y(0);
  NV_Ith_S(yvec, 1) = y(1);
  NV_Ith_S(yvec, 2) = y(2);
  const Eigen::Vector2d p = Eigen::Vector2d::Random(2);

  { // test evaluate
    const std::vector<boost::any>& result = root->Evaluate(yvec, p);
    
    EXPECT_DOUBLE_EQ(boost::any_cast<const double>(result[0]), p(0));
    EXPECT_DOUBLE_EQ(boost::any_cast<const double>(result[1]), y(2)*y(2)-p(1));
  }

  /*{ // test jacobian
    // input 0, output 0
    const boost::any jac00 = root->Jacobian(0, 0, y, p);
    const Eigen::MatrixXd& jac00ref = boost::any_cast<const Eigen::MatrixXd>(jac00);

    EXPECT_EQ(jac00ref.rows(), 1);
    EXPECT_EQ(jac00ref.cols(), 3);

    // input 0, output 1
    const boost::any jac01 = root->Jacobian(0, 1, y, p);
    const Eigen::MatrixXd& jac01ref = boost::any_cast<const Eigen::MatrixXd>(jac01);

    Eigen::MatrixXd expectedJac01 = Eigen::MatrixXd::Zero(1, 3);
    expectedJac01(0, 2) = 2.0*y(2);

    EXPECT_EQ(jac00ref.rows(), 1);
    EXPECT_EQ(jac00ref.cols(), 3);

    // check the values for input 0, output 0 and 1
    for( unsigned int j=0; j<3; ++j ) {
      EXPECT_DOUBLE_EQ(jac00ref(0,j), 0.0);

      EXPECT_DOUBLE_EQ(jac01ref(0,j), expectedJac01(0,j));
    }

    //  input 1, output 0
    const boost::any jac10 = root->Jacobian(1, 0, y, p);
    const Eigen::MatrixXd& jac10ref = boost::any_cast<const Eigen::MatrixXd>(jac10);

    Eigen::MatrixXd expectedJac10 = Eigen::MatrixXd::Zero(1, 2);
    expectedJac10(0, 0) = 1.0;

    EXPECT_EQ(jac10ref.rows(), 1);
    EXPECT_EQ(jac10ref.cols(), 2);

    // input 1, output 1
    const boost::any jac11 = root->Jacobian(1, 1, y, p);
    const Eigen::MatrixXd& jac11ref = boost::any_cast<const Eigen::MatrixXd>(jac11);

    Eigen::MatrixXd expectedJac11 = Eigen::MatrixXd::Zero(1, 2);
    expectedJac11(0, 1) = -1.0;

    EXPECT_EQ(jac11ref.rows(), 1);
    EXPECT_EQ(jac11ref.cols(), 2);

    // check the values of input 1, output 0 and 1
    for( unsigned int j=0; j<2; ++j ) {
      EXPECT_DOUBLE_EQ(jac10ref(0,j), expectedJac10(0,j));

      EXPECT_DOUBLE_EQ(jac11ref(0,j), expectedJac11(0,j));
    }
    }*/
}

TEST(RootfindingIVP, KineticsProblem) {
  // the right hand side of the ODE
  auto rhs = std::make_shared<KineticsRHS>();

  // integrate the ode until we find the root of this function
  auto root = std::make_shared<RootFunction>();

  /// Options for the ODE integrator
  pt::ptree pt;
  pt.put<double>("ODESolver.RelativeTolerance", 1.0e-10);
  pt.put<double>("ODESolver.AbsoluteTolerance", 1.0e-10);
  pt.put<double>("ODESolver.MaxStepSize", 1.0);
  pt.put<std::string>("ODESolver.MultistepMethod", "BDF");
  pt.put<std::string>("ODESolver.LinearSolver", "Dense");
  pt.put<std::string>("ODESolver.NonlinearSolver", "Newton");
  
  // the root finder
  auto rootfinder = std::make_shared<RootfindingIVP>(rhs, root, pt);

  // the input and output number is unknown
  EXPECT_EQ(rootfinder->numInputs, -1); // there is an optional input so even though the inputs to rhs and root are known
  EXPECT_EQ(rootfinder->numOutputs, -1); // there is an optional output that depends on the optional input

  // initial condition
  N_Vector ic = N_VNew_Serial(3);
  NV_Ith_S(ic, 0) = 1.0;
  NV_Ith_S(ic, 1) = 1.0;
  NV_Ith_S(ic, 2) = 1.0;

  // rhs parameters
  const Eigen::Vector2d a(1.0, 2.0);
  const double b = 1.0;

  // root parameters
  const Eigen::Vector2d p(5.0, 2.0);

  // evaluate the rootfinder
  const std::vector<boost::any>& result = rootfinder->Evaluate(ic, a, b, p);
  const N_Vector& rt = boost::any_cast<const N_Vector&>(result[0]);
  
  // the expected root
  const Eigen::VectorXd rtExpected = Eigen::Vector3d(1.5131802509043677, 2.7908899272796677, 1.4142135623730951);

  EXPECT_EQ(NV_LENGTH_S(rt), 3);
  for( unsigned int i=0; i<3; ++i ) {
    EXPECT_NEAR(NV_Ith_S(rt, i), rtExpected(i), 1.0e-8);
  }
}
