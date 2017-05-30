#include "gtest/gtest.h"

#include <Eigen/Core>

#include "MUQ/Modeling/ODE.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;

class RHS : public WorkPiece {
public:

  /// Constructor
  inline RHS() : WorkPiece(std::vector<std::string>({typeid(N_Vector).name(), typeid(double).name()}), std::vector<std::string>({typeid(N_Vector).name()})) {}

  inline virtual ~RHS() {}
  
private:

  inline virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override {
    // get the state vector
    const N_Vector& state = boost::any_cast<const N_Vector>(inputs[0]);

    // get the parameter (spring constant)
    const double k = boost::any_cast<const double>(inputs[1]);

    // set the output
    outputs.resize(1);
    outputs[0] = N_VNew_Serial(2);
    N_Vector& outref = boost::any_cast<N_Vector&>(outputs[0]);

    NV_Ith_S(outref, 0) = NV_Ith_S(state, 1);
    NV_Ith_S(outref, 1) = -k*NV_Ith_S(state, 0);
  }

  inline virtual void JacobianActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) override {
    // there is only one output
    assert(wrtOut==0);

    // get a reference to the vector
    const N_Vector& vecref = boost::any_cast<const N_Vector>(vec);

    if( wrtIn==0 ) { // wrt the state
      assert(NV_LENGTH_S(vecref)==2);
      
      // get the parameter (spring constant)
      const double k = boost::any_cast<const double>(inputs[1]);

      jacobianAction = N_VNew_Serial(2);
      N_Vector& jacAct = boost::any_cast<N_Vector&>(*jacobianAction);

      NV_Ith_S(jacAct, 0) = NV_Ith_S(vecref, 1);
      NV_Ith_S(jacAct, 1) = -k*NV_Ith_S(vecref, 0);
    } else if( wrtIn==1 ) {
      assert(NV_LENGTH_S(vecref)==1);
      
      // get the state vector
      const N_Vector& state = boost::any_cast<const N_Vector>(inputs[0]);

      jacobianAction = N_VNew_Serial(2);
      N_Vector& jacAct = boost::any_cast<N_Vector&>(*jacobianAction);

      NV_Ith_S(jacAct, 0) = 0.0;
      NV_Ith_S(jacAct, 1) = -NV_Ith_S(state, 0)*NV_Ith_S(vecref, 0);
    }
  }
  
  inline virtual void JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) override {
    // there is only one output
    assert(wrtOut==0);

    if( wrtIn==0 ) { // wrt the state
      // get the parameter (spring constant)
      const double k = boost::any_cast<const double>(inputs[1]);

      jacobian = NewDenseMat(2, 2);
      DlsMat& jac = boost::any_cast<DlsMat&>(*jacobian);
      SetToZero(jac);

      DENSE_ELEM(jac, 0, 1) = 1.0;
      DENSE_ELEM(jac, 1, 0) = -k;
    } else if( wrtIn==1 ) {
      // get the state vector
      const N_Vector& state = boost::any_cast<const N_Vector>(inputs[0]);

      jacobian = NewDenseMat(2, 1);
      DlsMat& jac = boost::any_cast<DlsMat&>(*jacobian);
      SetToZero(jac);

      DENSE_ELEM(jac, 1, 0) = -NV_Ith_S(state, 0);
    }
  }
};

// test the right hand side WorkPiece
TEST(ODEExample, RHS) {
  // create the right hand side
  auto rhs = std::make_shared<RHS>();

  // inputs
  N_Vector state = N_VNew_Serial(2);
  NV_Ith_S(state, 0) = 1.2;
  NV_Ith_S(state, 1) = 0.8;
  const double k = 2.25;

  { // test evaluate
    const std::vector<boost::any>& result = rhs->Evaluate(state, k);
    EXPECT_EQ(result.size(), 1);
    const N_Vector& rhsvec = boost::any_cast<const N_Vector>(result[0]);

    EXPECT_EQ(NV_LENGTH_S(rhsvec), 2);
    EXPECT_EQ(NV_Ith_S(rhsvec, 0), NV_Ith_S(state, 1));
    EXPECT_EQ(NV_Ith_S(rhsvec, 1), -k*NV_Ith_S(state, 0));
  }

  // expected jacobians
  Eigen::MatrixXd expectedJac0 = Eigen::MatrixXd::Zero(2, 2);
  expectedJac0(0, 1) = 1.0;
  expectedJac0(1, 0) = -k;

  Eigen::MatrixXd expectedJac1 = Eigen::MatrixXd::Zero(2, 1);
  expectedJac1(1, 0) = -NV_Ith_S(state, 0);
  
  { // test jacobian
    const boost::any& jac0 = rhs->Jacobian(0, 0, state, k);
    const DlsMat& jac0ref = boost::any_cast<const DlsMat>(jac0);

    EXPECT_EQ(jac0ref->M, 2); // check the number of rows
    EXPECT_EQ(jac0ref->N, 2); // check the number of cols
    for( unsigned int i=0; i<2; ++i ) {
      for( unsigned int j=0; j<2; ++j ) {
	EXPECT_DOUBLE_EQ(DENSE_ELEM(jac0ref, i, j), expectedJac0(i,j));
      }
    }
    
    const boost::any& jac1 = rhs->Jacobian(1, 0, state, k);
    const DlsMat& jac1ref = boost::any_cast<const DlsMat>(jac1);

    EXPECT_EQ(jac1ref->M, 2); // check the number of rows
    EXPECT_EQ(jac1ref->N, 1); // check the number of cols
    for( unsigned int i=0; i<2; ++i ) {
      EXPECT_DOUBLE_EQ(DENSE_ELEM(jac1ref, i, 0), expectedJac1(i,0));
    }
  }

  
  { // test jacobian action
    const Eigen::Vector2d vec0eigen = Eigen::VectorXd::Random(2);
    
    // a vector
    N_Vector vec0 = N_VNew_Serial(2);
    NV_Ith_S(vec0, 0) = vec0eigen(0);
    NV_Ith_S(vec0, 1) = vec0eigen(1);
    
    const boost::any& jacact0 = rhs->JacobianAction(0, 0, vec0, state, k);
    const N_Vector& jacact0ref = boost::any_cast<const N_Vector>(jacact0);

    // expected value
    const Eigen::Vector2d expectedJacact0 = expectedJac0*vec0eigen;

    EXPECT_EQ(NV_LENGTH_S(jacact0ref), 2); 
    for( unsigned int i=0; i<2; ++i ) {
      EXPECT_DOUBLE_EQ(NV_Ith_S(jacact0ref, i), expectedJacact0(i));
    }

    const Eigen::VectorXd vec1eigen = Eigen::VectorXd::Random(1);
    
    // a vector
    N_Vector vec1 = N_VNew_Serial(1);
    NV_Ith_S(vec1, 0) = vec1eigen(0);

    const boost::any& jacact1 = rhs->JacobianAction(1, 0, vec1, state, k);
    const N_Vector& jacact1ref = boost::any_cast<const N_Vector>(jacact1);

    // expected value
    const Eigen::Vector2d expectedJacact1 = expectedJac1*vec1eigen;

    EXPECT_EQ(NV_LENGTH_S(jacact1ref), 2); 
    for( unsigned int i=0; i<2; ++i ) {
      EXPECT_DOUBLE_EQ(NV_Ith_S(jacact1ref, i), expectedJacact1(i));
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

    // set the intial conditions
    ic = N_VNew_Serial(2);
    NV_Ith_S(ic, 0) = 1.0;
    NV_Ith_S(ic, 1) = 0.0;

    // solver options
    pt.put<double>("ODESolver.RelativeTolerance", 1.0e-8);
    pt.put<double>("ODESolver.AbsoluteTolerance", 1.0e-8);
    pt.put<double>("ODESolver.MaxStepSize", 1.0);
  }

  /// Default destructor
  virtual ~ODETests() {}

  virtual void TearDown() override {
    // the input/output number are unknown
    EXPECT_EQ(ode->numInputs, -1); // some of the inputs are the times where we need the state's value
    EXPECT_EQ(ode->numInputs, -1); // the outputs are the states at the specified times
    
    // integrate the ODE
    const std::vector<boost::any>& result = ode->Evaluate(ic, k, outTimes0);
    
    // check the result for the first vector of times
    const std::vector<N_Vector>& times0_state = boost::any_cast<const std::vector<N_Vector>&>(result[0]);
    EXPECT_EQ(times0_state.size(), outTimes0.size());
    for( unsigned int i=0; i<outTimes0.size(); ++i ) {
      EXPECT_NEAR(NV_Ith_S(times0_state[i], 0), cos(sqrt(k)*outTimes0(i)), 1.0e-6);
    }
  }

  /// The right hand side
  std::shared_ptr<RHS> rhs;

  /// The ode
  std::shared_ptr<ODE> ode;

  /// Options for the ODE integrator
  pt::ptree pt;

  /// The initial condition
  N_Vector ic;

  /// The spring constant
  const double k = 0.12;

  /// The output times 
  const Eigen::VectorXd outTimes0 = Eigen::VectorXd::LinSpaced(10, 0.0, 2.0);

private:
};

TEST_F(ODETests, BDFNewtonMethod) {
  pt.put<std::string>("ODESolver.MultistepMethod", "BDF");
  pt.put<std::string>("ODESolver.NonlinearSolver", "Newton");
  pt.put<std::string>("ODESolver.LinearSolver", "Dense");
  
  // create the ODE integrator
  ode = std::make_shared<ODE>(rhs, pt);

  // we need the state at these times
  const Eigen::Vector3d outTimes1(0.0, 0.5, 1.0);
  const double t0 = 1.0;
  const double t1 = 2.0;

  // integrate the ODE
  const std::vector<boost::any>& result = ode->Evaluate(ic, k, outTimes0, outTimes1, t0, t1);

  // check the result for the first vector of times
  const std::vector<N_Vector>& times0_state = boost::any_cast<const std::vector<N_Vector>&>(result[0]);
  EXPECT_EQ(times0_state.size(), outTimes0.size());
  for( unsigned int i=0; i<outTimes0.size(); ++i ) {
    EXPECT_NEAR(NV_Ith_S(times0_state[i], 0), cos(sqrt(k)*outTimes0(i)), 1.0e-6);
  }

  // check the result for the second vector of times
  const std::vector<N_Vector>& times1_state = boost::any_cast<const std::vector<N_Vector>&>(result[1]);
  EXPECT_EQ(times1_state.size(), outTimes1.size());
  for( unsigned int i=0; i<outTimes1.size(); ++i ) {
    EXPECT_NEAR(NV_Ith_S(times1_state[i], 0), cos(sqrt(k)*outTimes1(i)), 1.0e-6);
  }

  // check the result for the first scalar output
  const N_Vector& t0_state = boost::any_cast<const N_Vector&>(result[2]);
  EXPECT_NEAR(NV_Ith_S(t0_state, 0), cos(sqrt(k)*t0), 1.0e-6);

  // check the result for the second scalar output
  const N_Vector& t1_state = boost::any_cast<const N_Vector&>(result[3]);
  EXPECT_NEAR(NV_Ith_S(t1_state, 0), cos(sqrt(k)*t1), 1.0e-6);
}

TEST_F(ODETests, BDFIterMethod) {
  pt.put<std::string>("ODESolver.MultistepMethod", "BDF");
  pt.put<std::string>("ODESolver.NonlinearSolver", "Iter");
  pt.put<std::string>("ODESolver.LinearSolver", "Dense");
  
  // create the ODE integrator
  ode = std::make_shared<ODE>(rhs, pt);
}

TEST_F(ODETests, AdamsNewtonMethod) {
  pt.put<std::string>("ODESolver.MultistepMethod", "Adams");
  pt.put<std::string>("ODESolver.NonlinearSolver", "Newton");
  pt.put<std::string>("ODESolver.LinearSolver", "Dense");
  
  // create the ODE integrator
  ode = std::make_shared<ODE>(rhs, pt);
}

TEST_F(ODETests, AdamsIterMethod) {
  pt.put<std::string>("ODESolver.MultistepMethod", "Adams");
  pt.put<std::string>("ODESolver.NonlinearSolver", "Iter");
  pt.put<std::string>("ODESolver.LinearSolver", "Dense");
  
  // create the ODE integrator
  ode = std::make_shared<ODE>(rhs, pt);
}

TEST_F(ODETests, SPGMR) {
  pt.put<std::string>("ODESolver.MultistepMethod", "BDF");
  pt.put<std::string>("ODESolver.NonlinearSolver", "Newton");
  pt.put<std::string>("ODESolver.LinearSolver", "SPGMR");
  
  // create the ODE integrator
  ode = std::make_shared<ODE>(rhs, pt);
}

TEST_F(ODETests, SPBCG) {
  pt.put<std::string>("ODESolver.MultistepMethod", "Adams");
  pt.put<std::string>("ODESolver.NonlinearSolver", "Newton");
  pt.put<std::string>("ODESolver.LinearSolver", "SPBCG");
  
  // create the ODE integrator
  ode = std::make_shared<ODE>(rhs, pt);
}

TEST_F(ODETests, SPTFQMR) {
  pt.put<std::string>("ODESolver.MultistepMethod", "Adams");
  pt.put<std::string>("ODESolver.NonlinearSolver", "Iter");
  pt.put<std::string>("ODESolver.LinearSolver", "SPTFQMR");
  
  // create the ODE integrator
  ode = std::make_shared<ODE>(rhs, pt);
}
