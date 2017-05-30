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
    NV_Ith_S(ic, 0) = ic0;
    NV_Ith_S(ic, 1) = ic1;

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

    // forward sensitivity
    const boost::any& jac0 = ode->Jacobian(0, 0, ic, k, outTimes0);
    const std::vector<DlsMat>& jac0ref = boost::any_cast<const std::vector<DlsMat>&>(jac0);
    EXPECT_EQ(jac0ref.size(), outTimes0.size());
    
    const boost::any& jac1 = ode->Jacobian(1, 0, ic, k, outTimes0);
    const std::vector<DlsMat>& jac1ref = boost::any_cast<const std::vector<DlsMat>&>(jac1);
    EXPECT_EQ(jac1ref.size(), outTimes0.size());

    const boost::any& jac2 = ode->Jacobian(2, 0, ic, k, outTimes0);
    const std::vector<DlsMat>& jac2ref = boost::any_cast<const std::vector<DlsMat>&>(jac2);
    EXPECT_EQ(jac2ref.size(), outTimes0.size());

    for( unsigned int i=0; i<outTimes0.size(); ++i ) {
      const double time = outTimes0(i);

      // check evaluate values
      EXPECT_NEAR(NV_Ith_S(times0_state[i], 0), ic0*std::cos(std::sqrt(k)*time)+ic1/std::sqrt(k)*std::sin(std::sqrt(k)*time), 1.0e-6);
      EXPECT_NEAR(NV_Ith_S(times0_state[i], 1), -std::sqrt(k)*ic0*std::sin(std::sqrt(k)*time)+ic1*std::cos(std::sqrt(k)*time), 1.0e-6);

      // check jacobian wrt initial conditions
      EXPECT_EQ(jac0ref[i]->M, 2); // rows
      EXPECT_EQ(jac0ref[i]->N, 2); // cols
      EXPECT_NEAR(DENSE_ELEM(jac0ref[i], 0, 0), std::cos(std::sqrt(k)*time), 1.0e-6);
      EXPECT_NEAR(DENSE_ELEM(jac0ref[i], 0, 1), std::sin(std::sqrt(k)*time)/std::sqrt(k), 1.0e-6);
      EXPECT_NEAR(DENSE_ELEM(jac0ref[i], 1, 0), -std::sqrt(k)*std::sin(std::sqrt(k)*time), 1.0e-6);
      EXPECT_NEAR(DENSE_ELEM(jac0ref[i], 1, 1), std::cos(std::sqrt(k)*time), 1.0e-6);

      // check jacobian wrt spring constant
      EXPECT_EQ(jac1ref[i]->M, 2); // rows
      EXPECT_EQ(jac1ref[i]->N, 1); // cols
      EXPECT_NEAR(DENSE_ELEM(jac1ref[i], 0, 0), 0.5*(-ic0/std::sqrt(k)*time*std::sin(std::sqrt(k)*time)+ic1/k*time*std::cos(std::sqrt(k)*time)-ic1/std::pow(k, 1.5)*std::sin(std::sqrt(k)*time)), 1.0e-6);
      EXPECT_NEAR(DENSE_ELEM(jac1ref[i], 1, 0), 0.5*(-ic0*time*std::cos(std::sqrt(k)*time)-ic0/std::sqrt(k)*std::sin(std::sqrt(k)*time)-ic1/std::sqrt(k)*time*std::sin(std::sqrt(k)*time)), 1.0e-6);

      // check jacobian wrt output times --- which should just return the right hand side (derivative of state wrt time)
      EXPECT_EQ(jac2ref[i]->M, 2); // rows
      EXPECT_EQ(jac2ref[i]->N, 1); // cols
      EXPECT_DOUBLE_EQ(DENSE_ELEM(jac2ref[i], 0, 0), NV_Ith_S(times0_state[i], 1));
      EXPECT_DOUBLE_EQ(DENSE_ELEM(jac2ref[i], 1, 0), -k*NV_Ith_S(times0_state[i], 0));
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

  /// Initial condition for y0
  const double ic0 = 2.0;

  /// Initial condition for y1
  const double ic1 = 0.5;

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

  // compute jacobians of the first output
  const boost::any& jac00 = ode->Jacobian(0, 0, ic, k, outTimes0, outTimes1, t0, t1);
  const std::vector<DlsMat>& jac00ref = boost::any_cast<const std::vector<DlsMat>&>(jac00);
  EXPECT_EQ(jac00ref.size(), outTimes0.size());
  const boost::any& jac10 = ode->Jacobian(1, 0, ic, k, outTimes0, outTimes1, t0, t1);
  const std::vector<DlsMat>& jac10ref = boost::any_cast<const std::vector<DlsMat>&>(jac10);
  EXPECT_EQ(jac10ref.size(), outTimes0.size());

  // check the result for the first vector of times
  const std::vector<N_Vector>& times0_state = boost::any_cast<const std::vector<N_Vector>&>(result[0]);
  EXPECT_EQ(times0_state.size(), outTimes0.size());
  for( unsigned int i=0; i<outTimes0.size(); ++i ) {
    const double time = outTimes0(i);

    // check evaluate values
    EXPECT_NEAR(NV_Ith_S(times0_state[i], 0), ic0*std::cos(std::sqrt(k)*time)+ic1/std::sqrt(k)*std::sin(std::sqrt(k)*time), 1.0e-6);
    EXPECT_NEAR(NV_Ith_S(times0_state[i], 1), -std::sqrt(k)*ic0*std::sin(std::sqrt(k)*time)+ic1*std::cos(std::sqrt(k)*time), 1.0e-6);

    // check jacobian wrt initial conditions
    EXPECT_EQ(jac00ref[i]->M, 2); // rows
    EXPECT_EQ(jac00ref[i]->N, 2); // cols
    EXPECT_NEAR(DENSE_ELEM(jac00ref[i], 0, 0), std::cos(std::sqrt(k)*time), 1.0e-6);
    EXPECT_NEAR(DENSE_ELEM(jac00ref[i], 0, 1), std::sin(std::sqrt(k)*time)/std::sqrt(k), 1.0e-6);
    EXPECT_NEAR(DENSE_ELEM(jac00ref[i], 1, 0), -std::sqrt(k)*std::sin(std::sqrt(k)*time), 1.0e-6);
    EXPECT_NEAR(DENSE_ELEM(jac00ref[i], 1, 1), std::cos(std::sqrt(k)*time), 1.0e-6);
    
    // check jacobian wrt spring constant
    EXPECT_EQ(jac10ref[i]->M, 2); // rows
    EXPECT_EQ(jac10ref[i]->N, 1); // cols
    EXPECT_NEAR(DENSE_ELEM(jac10ref[i], 0, 0), 0.5*(-ic0/std::sqrt(k)*time*std::sin(std::sqrt(k)*time)+ic1/k*time*std::cos(std::sqrt(k)*time)-ic1/std::pow(k, 1.5)*std::sin(std::sqrt(k)*time)), 1.0e-6);
    EXPECT_NEAR(DENSE_ELEM(jac10ref[i], 1, 0), 0.5*(-ic0*time*std::cos(std::sqrt(k)*time)-ic0/std::sqrt(k)*std::sin(std::sqrt(k)*time)-ic1/std::sqrt(k)*time*std::sin(std::sqrt(k)*time)), 1.0e-6);
  }

  // compute jacobians of the second output
  const boost::any& jac01 = ode->Jacobian(0, 1, ic, k, outTimes0, outTimes1, t0, t1);
  const std::vector<DlsMat>& jac01ref = boost::any_cast<const std::vector<DlsMat>&>(jac01);
  EXPECT_EQ(jac01ref.size(), outTimes1.size());
  const boost::any& jac11 = ode->Jacobian(1, 1, ic, k, outTimes0, outTimes1, t0, t1);
  const std::vector<DlsMat>& jac11ref = boost::any_cast<const std::vector<DlsMat>&>(jac11);
  EXPECT_EQ(jac11ref.size(), outTimes1.size());

  // check the result for the second vector of times
  const std::vector<N_Vector>& times1_state = boost::any_cast<const std::vector<N_Vector>&>(result[1]);
  EXPECT_EQ(times1_state.size(), outTimes1.size());
  for( unsigned int i=0; i<outTimes1.size(); ++i ) {
    const double time = outTimes1(i);
    
    EXPECT_NEAR(NV_Ith_S(times1_state[i], 0), ic0*std::cos(std::sqrt(k)*time)+ic1/std::sqrt(k)*std::sin(std::sqrt(k)*time), 1.0e-6);
    EXPECT_NEAR(NV_Ith_S(times1_state[i], 1), -std::sqrt(k)*ic0*std::sin(std::sqrt(k)*time)+ic1*std::cos(std::sqrt(k)*time), 1.0e-6);

    // check jacobian wrt initial conditions
    EXPECT_EQ(jac01ref[i]->M, 2); // rows
    EXPECT_EQ(jac01ref[i]->N, 2); // cols
    EXPECT_NEAR(DENSE_ELEM(jac01ref[i], 0, 0), std::cos(std::sqrt(k)*time), 1.0e-6);
    EXPECT_NEAR(DENSE_ELEM(jac01ref[i], 0, 1), std::sin(std::sqrt(k)*time)/std::sqrt(k), 1.0e-6);
    EXPECT_NEAR(DENSE_ELEM(jac01ref[i], 1, 0), -std::sqrt(k)*std::sin(std::sqrt(k)*time), 1.0e-6);
    EXPECT_NEAR(DENSE_ELEM(jac01ref[i], 1, 1), std::cos(std::sqrt(k)*time), 1.0e-6);
    
    // check jacobian wrt spring constant
    EXPECT_EQ(jac11ref[i]->M, 2); // rows
    EXPECT_EQ(jac11ref[i]->N, 1); // cols
    EXPECT_NEAR(DENSE_ELEM(jac11ref[i], 0, 0), 0.5*(-ic0/std::sqrt(k)*time*std::sin(std::sqrt(k)*time)+ic1/k*time*std::cos(std::sqrt(k)*time)-ic1/std::pow(k, 1.5)*std::sin(std::sqrt(k)*time)), 1.0e-6);
    EXPECT_NEAR(DENSE_ELEM(jac11ref[i], 1, 0), 0.5*(-ic0*time*std::cos(std::sqrt(k)*time)-ic0/std::sqrt(k)*std::sin(std::sqrt(k)*time)-ic1/std::sqrt(k)*time*std::sin(std::sqrt(k)*time)), 1.0e-6);
  }
  
  // check the result for the first scalar output
  const N_Vector& t0_state = boost::any_cast<const N_Vector&>(result[2]);
  EXPECT_NEAR(NV_Ith_S(t0_state, 0), ic0*std::cos(std::sqrt(k)*t0)+ic1/std::sqrt(k)*std::sin(std::sqrt(k)*t0), 1.0e-6);
  EXPECT_NEAR(NV_Ith_S(t0_state, 1), -std::sqrt(k)*ic0*std::sin(std::sqrt(k)*t0)+ic1*std::cos(std::sqrt(k)*t0), 1.0e-6);

  // compute jacobians of the first scalar output
  const boost::any& jac02 = ode->Jacobian(0, 2, ic, k, outTimes0, outTimes1, t0, t1);
  const DlsMat& jac02ref = boost::any_cast<const DlsMat&>(jac02);

  EXPECT_EQ(jac02ref->M, 2); // rows
  EXPECT_EQ(jac02ref->N, 2); // cols
  EXPECT_NEAR(DENSE_ELEM(jac02ref, 0, 0), std::cos(std::sqrt(k)*t0), 1.0e-6);
  EXPECT_NEAR(DENSE_ELEM(jac02ref, 0, 1), std::sin(std::sqrt(k)*t0)/std::sqrt(k), 1.0e-6);
  EXPECT_NEAR(DENSE_ELEM(jac02ref, 1, 0), -std::sqrt(k)*std::sin(std::sqrt(k)*t0), 1.0e-6);
  EXPECT_NEAR(DENSE_ELEM(jac02ref, 1, 1), std::cos(std::sqrt(k)*t0), 1.0e-6);

  const boost::any& jac12 = ode->Jacobian(1, 2, ic, k, outTimes0, outTimes1, t0, t1);
  const DlsMat& jac12ref = boost::any_cast<const DlsMat&>(jac12);

  EXPECT_EQ(jac12ref->M, 2); // rows
  EXPECT_EQ(jac12ref->N, 1); // cols
  EXPECT_NEAR(DENSE_ELEM(jac12ref, 0, 0), 0.5*(-ic0/std::sqrt(k)*t0*std::sin(std::sqrt(k)*t0)+ic1/k*t0*std::cos(std::sqrt(k)*t0)-ic1/std::pow(k, 1.5)*std::sin(std::sqrt(k)*t0)), 1.0e-6);
  EXPECT_NEAR(DENSE_ELEM(jac12ref, 1, 0), 0.5*(-ic0*t0*std::cos(std::sqrt(k)*t0)-ic0/std::sqrt(k)*std::sin(std::sqrt(k)*t0)-ic1/std::sqrt(k)*t0*std::sin(std::sqrt(k)*t0)), 1.0e-6);

  // check the result for the second scalar output
  const N_Vector& t1_state = boost::any_cast<const N_Vector&>(result[3]);
  EXPECT_NEAR(NV_Ith_S(t1_state, 0), ic0*std::cos(std::sqrt(k)*t1)+ic1/std::sqrt(k)*std::sin(std::sqrt(k)*t1), 1.0e-6);
  EXPECT_NEAR(NV_Ith_S(t1_state, 1), -std::sqrt(k)*ic0*std::sin(std::sqrt(k)*t1)+ic1*std::cos(std::sqrt(k)*t1), 1.0e-6);

  // compute jacobians of the second scalar output
  const boost::any& jac03 = ode->Jacobian(0, 3, ic, k, outTimes0, outTimes1, t0, t1);
  const DlsMat& jac03ref = boost::any_cast<const DlsMat&>(jac03);

  EXPECT_EQ(jac03ref->M, 2); // rows
  EXPECT_EQ(jac03ref->N, 2); // cols
  EXPECT_NEAR(DENSE_ELEM(jac03ref, 0, 0), std::cos(std::sqrt(k)*t1), 1.0e-6);
  EXPECT_NEAR(DENSE_ELEM(jac03ref, 0, 1), std::sin(std::sqrt(k)*t1)/std::sqrt(k), 1.0e-6);
  EXPECT_NEAR(DENSE_ELEM(jac03ref, 1, 0), -std::sqrt(k)*std::sin(std::sqrt(k)*t1), 1.0e-6);
  EXPECT_NEAR(DENSE_ELEM(jac03ref, 1, 1), std::cos(std::sqrt(k)*t1), 1.0e-6);

  const boost::any& jac13 = ode->Jacobian(1, 3, ic, k, outTimes0, outTimes1, t0, t1);
  const DlsMat& jac13ref = boost::any_cast<const DlsMat&>(jac13);

  EXPECT_EQ(jac13ref->M, 2); // rows
  EXPECT_EQ(jac13ref->N, 1); // cols
  EXPECT_NEAR(DENSE_ELEM(jac13ref, 0, 0), 0.5*(-ic0/std::sqrt(k)*t1*std::sin(std::sqrt(k)*t1)+ic1/k*t1*std::cos(std::sqrt(k)*t1)-ic1/std::pow(k, 1.5)*std::sin(std::sqrt(k)*t1)), 1.0e-6);
  EXPECT_NEAR(DENSE_ELEM(jac13ref, 1, 0), 0.5*(-ic0*t1*std::cos(std::sqrt(k)*t1)-ic0/std::sqrt(k)*std::sin(std::sqrt(k)*t1)-ic1/std::sqrt(k)*t1*std::sin(std::sqrt(k)*t1)), 1.0e-6);

  // check jacobian wrt output times (scalar case) --- which should just return the right hand side (derivative of state wrt time)
  const boost::any& jac23 = ode->Jacobian(2, 3, ic, k, outTimes0, outTimes1, t0, t1);
  const DlsMat& jac23ref = boost::any_cast<const DlsMat&>(jac23);

  EXPECT_EQ(jac23ref->M, 2); // rows
  EXPECT_EQ(jac23ref->N, 1); // cols
  EXPECT_NEAR(DENSE_ELEM(jac23ref, 0, 0), NV_Ith_S(t1_state, 1), 1.0e-6);
  EXPECT_NEAR(DENSE_ELEM(jac23ref, 1, 0), -k*NV_Ith_S(t1_state, 0), 1.0e-6);
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
