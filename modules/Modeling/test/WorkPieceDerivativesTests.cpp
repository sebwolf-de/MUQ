#include "gtest/gtest.h"

#include <Eigen/Core>

#include "MUQ/Modeling/WorkPiece.h"

using namespace muq::Modeling;

/// A linear function to test the finite difference derivatives
class Linear : public WorkPiece {
public:

  Linear(Eigen::MatrixXd const& Q, Eigen::VectorXd const& b) :
    WorkPiece(std::vector<std::string>({typeid(double).name(), typeid(Eigen::VectorXd).name()}), std::vector<std::string>({typeid(std::string).name(), typeid(Eigen::VectorXd).name()})),
    Q(Q), b(b)
  {}

  virtual ~Linear() {}

private:

  virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override {
    outputs.resize(2);

    // a double to scale the linear part
    const double a = boost::any_cast<double>(inputs[0]);

    // constant reference to the input vector
    const Eigen::VectorXd& in = boost::any_cast<const Eigen::VectorXd&>(inputs[1]);

    // the first output is a string
    outputs[0] = (std::string)"string";

    // compute the linear function
    outputs[1] = (Eigen::VectorXd)(Q*in + b);
  }

  /// Matrix for the linear part
  const Eigen::MatrixXd Q;

  /// Vector for the constant part
  const Eigen::VectorXd b;
};

/// A quadratic function to test the derivatives
class Quadratic : public WorkPiece {
public:

  Quadratic(Eigen::MatrixXd const& Q, Eigen::VectorXd const& a, Eigen::VectorXd const& b) :
    WorkPiece(std::vector<std::string>(1, typeid(Eigen::VectorXd).name()), std::vector<std::string>(1, typeid(double).name())),
    Q(Q), a(a), b(b)
  {}

  virtual ~Quadratic() {}
  
private:

  virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override {
    outputs.resize(1);

    // constant reference to the input vector
    const Eigen::VectorXd& in = boost::any_cast<Eigen::VectorXd>(inputs[0]);
    // compute the quadratic function
    outputs[0] = (in.transpose()*Q*in + a.transpose()*in + b) (0);
  }

  virtual void JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) override {
    // constant reference to the input vector
    const Eigen::VectorXd& in = boost::any_cast<const Eigen::VectorXd&>(inputs[0]);

    // compute the Jacobian
    jacobian = (Eigen::MatrixXd)(2.0*in.transpose()*Q + a.transpose());
  }

  virtual void JacobianActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) override {
    // constant reference to the input vector
    const Eigen::VectorXd& in = boost::any_cast<const Eigen::VectorXd&>(inputs[0]);

    // constant reference to the vector we are applying the Jacobian to 
    const Eigen::VectorXd& appvec = boost::any_cast<const Eigen::VectorXd&>(vec);

    // compute the Jacobian
    jacobianAction = (2.0*in.transpose()*Q*appvec + a.transpose()*appvec) (0);
  }

  virtual void JacobianTransposeActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) override {
    // constant reference to the input vector
    const Eigen::VectorXd& in = boost::any_cast<const Eigen::VectorXd&>(inputs[0]);

    // constant reference to the vector we are applying the Jacobian to 
    const Eigen::VectorXd& appvec = boost::any_cast<const Eigen::VectorXd&>(vec);

    // compute the Jacobian
    jacobianTransposeAction = (Eigen::VectorXd)(2.0*Q.transpose()*in*appvec + a*appvec);
  }

  /// Matrix for the quadratic part
  const Eigen::MatrixXd Q;

  /// Vector for the linear part
  const Eigen::VectorXd a;

  /// Vector for the constant part
  const Eigen::VectorXd b;
};

/// A class to test the behavior of WorkPiece with various input/output types/numbers
class WorkPieceDerivativesTests : public::testing::Test {
public:

  /// Default constructor
  WorkPieceDerivativesTests() {    
    // a random matrix (for the quadratic term)
    Q = Eigen::MatrixXd::Random(N, N);

    // a random vector (for the linear term)
    a = Eigen::VectorXd::Random(N);

    // a random vector (for the constant term)
    b = Eigen::VectorXd::Random(1);

    // create a linear polynomial
    lin = std::make_shared<Linear>(Q, a);

    // create a quadratic polynomial
    quad = std::make_shared<Quadratic>(Q, a, b);
  }

  /// Default destructor
  virtual ~WorkPieceDerivativesTests() {}

  /// The size of the system
  const unsigned int N = 5;

  /// Matrix for the quadratic part
  Eigen::MatrixXd Q;

  /// Vector for the linear part
  Eigen::VectorXd a;

  /// Vector for the constant part
  Eigen::VectorXd b;  

  /// A quadratic function
  std::shared_ptr<Quadratic> quad;

  /// A linear function
  std::shared_ptr<Linear> lin;

private:
};

TEST_F(WorkPieceDerivativesTests, LinearFunction) {
  // check input/output sizes
  EXPECT_EQ(lin->numInputs, 2);
  EXPECT_EQ(lin->numOutputs, 2);

  // choose random inputs
  const Eigen::VectorXd in = Eigen::VectorXd::Random(N);
  const double scalar = 3.5;

  // evaluate 
  auto result = lin->Evaluate(scalar, in);

  // make sure we get the expected result
  EXPECT_EQ(result.size(), 2);
  EXPECT_EQ(boost::any_cast<std::string>(result[0]).compare("string"), 0);
  const Eigen::VectorXd& vec = boost::any_cast<Eigen::VectorXd>(result[1]);

  // the expected vector
  const Eigen::VectorXd expectedVec = Q*in+a;
  
  EXPECT_EQ(vec.size(), N);
  EXPECT_EQ(expectedVec.size(), N);
  for( unsigned int i=0; i<N; ++i ) {
    EXPECT_DOUBLE_EQ(vec(i), expectedVec(i));
  }

  // compute the Jacobian and get a reference to it
  auto jac = lin->Jacobian(1, 1, scalar, in);
  const Eigen::MatrixXd& jacref = boost::any_cast<const Eigen::MatrixXd&>(jac);

  EXPECT_EQ(jacref.rows(), N);
  EXPECT_EQ(jacref.cols(), N);
  for( unsigned int i=0; i<N; ++i ) {
    for( unsigned int j=0; j<N; ++j ) {
      // its linear so FD should be exact, but the error is very small ...
      EXPECT_NEAR(jacref(i,j), Q(i,j), 1.0e-9);
    }
  }
}

TEST_F(WorkPieceDerivativesTests, QuadraticFunction) {
  // check input/output sizes
  EXPECT_EQ(quad->numInputs, 1);
  EXPECT_EQ(quad->numOutputs, 1);

  // choose random inputs
  const Eigen::VectorXd in = Eigen::VectorXd::Random(N);

  // evaluate and make sure we get the expected result
  auto result = quad->Evaluate(in);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(result[0]), (in.transpose()*Q*in + a.transpose()*in + b) (0));

  // evaluate the jacobian 
  auto jacBoost = quad->Jacobian(0, 0, in);
  const Eigen::MatrixXd& jac = boost::any_cast<Eigen::MatrixXd>(jacBoost);

  // compute the expected jacobian
  const Eigen::MatrixXd jacExpected = 2.0*in.transpose()*Q + a.transpose();

  // make sure the jacobians match
  EXPECT_EQ(jac.rows(), 1);
  EXPECT_EQ(jacExpected.rows(), 1);
  EXPECT_EQ(jac.cols(), N);
  EXPECT_EQ(jacExpected.cols(), N);
  for( unsigned int i=0; i<N; ++i ) {
    EXPECT_DOUBLE_EQ(jac(0,i), jacExpected(0,i));
  }

  // a random vector to apply the Jacobian to
  const Eigen::VectorXd vec0 = Eigen::VectorXd::Random(N);

  // evaluate the jacobian action
  auto jacAction = quad->JacobianAction(0, 0, vec0, in);

  // compute the expected jacobian action
  const double jacActionExpected = (2.0*in.transpose()*Q*vec0 + a.transpose()*vec0) (0);

  // make sure the jacobian action matches
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(jacActionExpected), jacActionExpected);

  // a random vector to apply the Jacobian transpose to
  const Eigen::VectorXd vec1 = Eigen::VectorXd::Random(1);

  // evaluate the jacobian action
  auto jacTransposeActionBoost = quad->JacobianTransposeAction(0, 0, vec1, in);
  const Eigen::VectorXd& jacTransposeAction = boost::any_cast<Eigen::VectorXd>(jacTransposeActionBoost);

  // compute the expected jacobian action
  const Eigen::VectorXd jacTransposeActionExpected = 2.0*Q.transpose()*in*vec1 + a*vec1;

  // make sure the jacobian action matches
  EXPECT_EQ(jacTransposeAction.rows(), N);
  EXPECT_EQ(jacTransposeActionExpected.rows(), N);
  EXPECT_EQ(jacTransposeAction.cols(), 1);
  EXPECT_EQ(jacTransposeActionExpected.cols(), 1);
  for( unsigned int i=0; i<N; ++i ) {
    EXPECT_DOUBLE_EQ(jacTransposeAction(i,0), jacTransposeActionExpected(i,0));
  }
}
