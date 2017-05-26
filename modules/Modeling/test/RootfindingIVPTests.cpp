#include "gtest/gtest.h"

#include <Eigen/Core>

#include "MUQ/Modeling/RootfindingIVP.h"

using namespace muq::Modeling;

/// The right hand side for a chemical kenetics ODE (from Sundials example)
class KineticsRHS : public WorkPiece {
public:

  /// Constructor
  inline KineticsRHS() : WorkPiece(std::vector<std::string>({typeid(Eigen::Vector3d).name(), typeid(Eigen::Vector2d).name(), typeid(double).name()}), std::vector<std::string>(1, typeid(Eigen::Vector3d).name())) {}

  inline virtual ~KineticsRHS() {}
  
private:

  inline virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override {
    outputs.resize(1);
    
    // dependent variables
    const Eigen::Vector3d& y = boost::any_cast<const Eigen::Vector3d>(inputs[0]);

    // parameters
    const Eigen::Vector2d& a = boost::any_cast<const Eigen::Vector2d>(inputs[1]);
    const double b = boost::any_cast<const double>(inputs[2]);

    // derivative wrt time
    outputs[0] = Eigen::Vector3d(3);
    Eigen::Vector3d& refout = boost::any_cast<Eigen::Vector3d&>(outputs[0]);
    refout(0) = a(0)*y(0); // dy1dt
    refout(1) = a(1)*y(0)*y(1); // dy2dt
    refout(2) = b*b; // dy3dt
  }

  inline virtual void JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) override {
    // dependent variables
    const Eigen::Vector3d& y = boost::any_cast<const Eigen::Vector3d>(inputs[0]);

    // parameters
    const Eigen::Vector2d& a = boost::any_cast<const Eigen::Vector2d>(inputs[1]);
    const double b = boost::any_cast<const double>(inputs[2]);

    assert(wrtOut==0);
    if( wrtIn==0 ) { // dervative wrt y
      // get a reference to the jacobian matrix
      jacobian = (Eigen::MatrixXd)Eigen::MatrixXd::Zero(3, 3);
      Eigen::MatrixXd& jac = boost::any_cast<Eigen::MatrixXd&>(*jacobian);

      // implement the jacobian
      jac(0, 0) = a(0);
      jac(1, 0) = a(1)*y(1);
      jac(1, 1) = a(1)*y(0);
    } else if( wrtIn==1 ) { // derivative wrt a
      // get a reference to the jacobian matrix
      jacobian = (Eigen::MatrixXd)Eigen::MatrixXd::Zero(3, 2);
      Eigen::MatrixXd& jac = boost::any_cast<Eigen::MatrixXd&>(*jacobian);

      jac(0, 0) = y(0);
      jac(1, 1) = y(0)*y(1);
    } else if( wrtIn==2 ) { // dervative wrt b
      // get a reference to the jacobian matrix
      jacobian = (Eigen::MatrixXd)Eigen::MatrixXd::Zero(3, 1);
      Eigen::MatrixXd& jac = boost::any_cast<Eigen::MatrixXd&>(*jacobian);

      jac(2, 0) = 2.0*b;
    }
  }
};

/// We want to the integrate the kinetics probem until we find the root of this function (from a Sundials example)
class RootFunction : public WorkPiece {
public:

  /// Constructor
  inline RootFunction() : WorkPiece(std::vector<std::string>({typeid(Eigen::Vector3d).name(), typeid(Eigen::Vector2d).name()}), std::vector<std::string>(2, typeid(double).name())) {}

  inline virtual ~RootFunction() {}
  
private:

  inline virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override {
    // dependent variables
    const Eigen::Vector3d& y = boost::any_cast<const Eigen::Vector3d>(inputs[0]);

    // parameters
    const Eigen::Vector2d& p = boost::any_cast<const Eigen::Vector2d>(inputs[1]);

    // compute the root function
    outputs.resize(2);
    outputs[0] = p(0); // root 1 (this will only be zero if p(0) is ...)
    outputs[1] = y(2)*y(2)-p(1); // root 2 
  }

  inline virtual void JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) override {
    // dependent variables
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
    }
  }
};

// test the RHS of the kinetics problem (not actually part of MUQ, but a good sanity check)
TEST(KineticsProblemTest, RHS) {
  // the right hand side of the ODE
  auto rhs = std::make_shared<KineticsRHS>();

  // inputs
  const Eigen::Vector3d y = Eigen::Vector3d::Random(3);
  const Eigen::Vector2d a = Eigen::Vector2d::Random(2);
  const double b = 3.5;

  { // test evaluate
    const std::vector<boost::any>& result = rhs->Evaluate(y, a, b);
    const Eigen::Vector3d& rhs = boost::any_cast<const Eigen::Vector3d>(result[0]);

    // check result
    EXPECT_DOUBLE_EQ(rhs(0), a(0)*y(0));
    EXPECT_DOUBLE_EQ(rhs(1), a(1)*y(0)*y(1));
    EXPECT_DOUBLE_EQ(rhs(2), b*b);
  }

  { // test jacobian
    // wrt to the state y
    const boost::any& jac0 = rhs->Jacobian(0, 0, y, a, b);
    const Eigen::MatrixXd jac0ref = boost::any_cast<const Eigen::MatrixXd>(jac0);

    Eigen::MatrixXd expectedJac0 = Eigen::MatrixXd::Zero(3, 3);
    expectedJac0(0, 0) = a(0);
    expectedJac0(1, 0) = a(1)*y(1);
    expectedJac0(1, 1) = a(1)*y(0);

    EXPECT_EQ(jac0ref.rows(), 3);
    EXPECT_EQ(jac0ref.cols(), 3);

    // wrt to the parameters a
    const boost::any& jac1 = rhs->Jacobian(1, 0, y, a, b);
    const Eigen::MatrixXd jac1ref = boost::any_cast<const Eigen::MatrixXd>(jac1);

    Eigen::MatrixXd expectedJac1 = Eigen::MatrixXd::Zero(3, 2);
    expectedJac1(0, 0) = y(0);
    expectedJac1(1, 1) = y(0)*y(1);

    EXPECT_EQ(jac1ref.rows(), 3);
    EXPECT_EQ(jac1ref.cols(), 2);

    // wrt to the parameters b
    const boost::any& jac2 = rhs->Jacobian(2, 0, y, a, b);
    const Eigen::MatrixXd jac2ref = boost::any_cast<const Eigen::MatrixXd>(jac2);

    Eigen::MatrixXd expectedJac2 = Eigen::MatrixXd::Zero(3, 1);
    expectedJac2(2, 0) = 2.0*b;

    EXPECT_EQ(jac2ref.rows(), 3);
    EXPECT_EQ(jac2ref.cols(), 1);

    // check the values
    for( unsigned int i=0; i<3; ++i ) {
      for( unsigned int j=0; j<3; ++j ) {
	EXPECT_DOUBLE_EQ(jac0ref(i,j), expectedJac0(i,j));
      }

      for( unsigned int j=0; j<2; ++j ) {
	EXPECT_DOUBLE_EQ(jac1ref(i,j), expectedJac1(i,j));
      }

      EXPECT_DOUBLE_EQ(jac2ref(i,0), expectedJac2(i,0));
    }
  }
}

// test the root function of the kinetics problem (not actually part of MUQ, but a good sanity check)
TEST(KineticsProblemTest, Root) {
  // we are trying to find the root of this function
  auto root = std::make_shared<RootFunction>();

  // inputs
  const Eigen::Vector3d y = Eigen::Vector3d::Random(3);
  const Eigen::Vector2d p = Eigen::Vector2d::Random(2);

  { // test evaluate
    const std::vector<boost::any>& result = root->Evaluate(y, p);
    
    EXPECT_DOUBLE_EQ(boost::any_cast<const double>(result[0]), p(0));
    EXPECT_DOUBLE_EQ(boost::any_cast<const double>(result[1]), y(2)*y(2)-p(1));
  }

  { // test jacobian
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
  }
}

TEST(RootfindingIVP, KineticsProblem) {
  // the right hand side of the ODE
  auto rhs = std::make_shared<KineticsRHS>();

  // integrate the ode until we find the root of this function
  auto root = std::make_shared<RootFunction>();

  // the root finder
  auto rootfinder = std::make_shared<RootfindingIVP>(rhs, root);

  // the input and output number is unknown
  EXPECT_EQ(rootfinder->numInputs, -1); // there is an optional input so even though the inputs to rhs and root are known, this is -1
  EXPECT_EQ(rootfinder->numOutputs, -1); // there is an optional output that depends on the optional input

  // initial condition
  const Eigen::Vector3d ic(1.0, 1.0, 1.0);

  // rhs parameters
  const Eigen::Vector2d a(1.0, 2.0);
  const double b = 1.0;

  // root parameters
  const Eigen::Vector2d p(5.0, 2.0);

  // evaluate the rootfinder
  const std::vector<boost::any>& result = rootfinder->Evaluate(ic, a, b, p);
}
