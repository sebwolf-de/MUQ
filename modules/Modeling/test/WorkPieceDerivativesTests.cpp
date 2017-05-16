#include "gtest/gtest.h"

#include <Eigen/Core>

#include "MUQ/Modeling/WorkPiece.h"

using namespace muq::Modeling;

/// A polynomial function to test the derivatives
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
    const Eigen::VectorXd& in = boost::any_cast<Eigen::VectorXd>(inputs[0]);

    // compute the Jacobian
    jacobian = (Eigen::MatrixXd)(2.0*in.transpose()*Q + a.transpose());
  }

  virtual void JacobianActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) override {
    // constant reference to the input vector
    const Eigen::VectorXd& in = boost::any_cast<Eigen::VectorXd>(inputs[0]);

    // constant reference to the vector we are applying the Jacobian to 
    const Eigen::VectorXd& appvec = boost::any_cast<Eigen::VectorXd>(vec);

    // compute the Jacobian
    jacobianAction = (2.0*in.transpose()*Q*appvec + a.transpose()*appvec) (0);
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

    // create a quadratic polynomial
    poly = std::make_shared<Quadratic>(Q, a, b);
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

  /// A polynomial function
  std::shared_ptr<Quadratic> poly;
private:
};

TEST_F(WorkPieceDerivativesTests, QuadraticFunction) {
  // check input/output sizes
  EXPECT_EQ(poly->numInputs, 1);
  EXPECT_EQ(poly->numOutputs, 1);

  // choose random inputs
  const Eigen::VectorXd in = Eigen::VectorXd::Random(N);

  // evaluate and make sure we get the expected result
  auto result = poly->Evaluate(in);
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(result[0]), (in.transpose()*Q*in + a.transpose()*in + b) (0));

  // evaluate the jacobian 
  auto jacBoost = poly->Jacobian(0, 0, in);
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
  const Eigen::VectorXd vec = Eigen::VectorXd::Random(N);

  // evaluate the jacobian action
  auto jacAction = poly->JacobianAction(0, 0, vec, in);

  // compute the expected jacobian action
  const double jacActionExpected = (2.0*in.transpose()*Q*vec + a.transpose()*vec) (0);

  // make sure the jacobian action matches
  EXPECT_DOUBLE_EQ(boost::any_cast<double>(jacActionExpected), jacActionExpected);

}
