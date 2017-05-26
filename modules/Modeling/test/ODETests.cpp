#include "gtest/gtest.h"

#include <Eigen/Core>

#include "MUQ/Modeling/ODE.h"

using namespace muq::Modeling;

class RHS : public WorkPiece {
public:

  /// Constructor
  inline RHS() : WorkPiece(std::vector<std::string>({typeid(Eigen::Vector2d).name(), typeid(double).name()}), std::vector<std::string>({typeid(Eigen::Vector2d).name()})) {}

  inline virtual ~RHS() {}
  
private:

  inline virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override {
    // get the state vector
    const Eigen::VectorXd& state = boost::any_cast<const Eigen::Vector2d>(inputs[0]);

    // get the parameter (spring constant)
    const double k = boost::any_cast<const double>(inputs[1]);

    // set the output
    outputs.resize(1);
    outputs[0] = Eigen::Vector2d(0.0, 0.0);
    Eigen::Vector2d& outref = boost::any_cast<Eigen::Vector2d&>(outputs[0]);

    outref[0] = state(1);
    outref[1] = -k*state(0);
  }

  inline virtual void JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) override {
    // there is only one output
    assert(wrtOut==0);

    if( wrtIn==0 ) { // wrt the state
      // get the parameter (spring constant)
      const double k = boost::any_cast<const double>(inputs[1]);

      jacobian = (Eigen::MatrixXd)Eigen::MatrixXd::Zero(2, 2);
      Eigen::MatrixXd& jac = boost::any_cast<Eigen::MatrixXd&>(*jacobian);

      jac(0, 1) = 1.0;
      jac(1, 0) = -k;
    } else if( wrtIn==1 ) {
      // get the state vector
      const Eigen::VectorXd& state = boost::any_cast<const Eigen::Vector2d>(inputs[0]);

      jacobian = (Eigen::MatrixXd)Eigen::MatrixXd::Zero(2, 1);
      Eigen::MatrixXd& jac = boost::any_cast<Eigen::MatrixXd&>(*jacobian);
      
      jac(1, 0) = -state(0);
    }
  }
};

// test the right hand side WorkPiece
TEST(ODEExample, RHS) {
  // create the right hand side
  auto rhs = std::make_shared<RHS>();

  // inputs
  const Eigen::Vector2d state = Eigen::VectorXd::Random(2);
  const double k = 2.25;

  { // test evaluate
    const std::vector<boost::any>& result = rhs->Evaluate(state, k);
    EXPECT_EQ(result.size(), 1);
    const Eigen::Vector2d& rhsvec = boost::any_cast<const Eigen::Vector2d>(result[0]);

    EXPECT_EQ(rhsvec.size(), 2);
    EXPECT_EQ(rhsvec(0), state(1));
    EXPECT_EQ(rhsvec(1), -k*state(0));
  }

  { // test jacobian
    const boost::any& jac0 = rhs->Jacobian(0, 0, state, k);
    const Eigen::MatrixXd& jac0ref = boost::any_cast<const Eigen::MatrixXd>(jac0);

    Eigen::MatrixXd expectedJac0 = Eigen::MatrixXd::Zero(2, 2);
    expectedJac0(0, 1) = 1.0;
    expectedJac0(1, 0) = -k;

    EXPECT_EQ(jac0ref.rows(), 2);
    EXPECT_EQ(jac0ref.cols(), 2);
    for( unsigned int i=0; i<2; ++i ) {
      for( unsigned int j=0; j<2; ++j ) {
	EXPECT_DOUBLE_EQ(jac0ref(i,j), expectedJac0(i,j));
      }
    }
    
    const boost::any& jac1 = rhs->Jacobian(1, 0, state, k);
    const Eigen::MatrixXd& jac1ref = boost::any_cast<const Eigen::MatrixXd>(jac1);

    Eigen::MatrixXd expectedJac1 = Eigen::MatrixXd::Zero(2, 1);
    expectedJac1(1, 0) = -state(0);

    EXPECT_EQ(jac1ref.rows(), 2);
    EXPECT_EQ(jac1ref.cols(), 1);
    for( unsigned int i=0; i<2; ++i ) {
      EXPECT_DOUBLE_EQ(jac1ref(i,0), expectedJac1(i,0));
    }
  }
}

/// A class to test the behavior of WorkPiece with various input/output types/numbers
class ODETests : public::testing::Test {
public:

  /// Default constructor
  ODETests() {
    // create the right hand side
    rhs = std::make_shared<RHS>();
  }

  /// Default destructor
  virtual ~ODETests() {}

  std::shared_ptr<RHS> rhs;

private:
};

TEST_F(ODETests, DenseSolver) {
  // create the ODE integrator
  auto ode = std::make_shared<ODE>(rhs);
}

