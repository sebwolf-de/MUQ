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

  virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override {
    outputs.resize(1);

    // constant reference to the input vector
    const Eigen::VectorXd& in = boost::any_cast<Eigen::VectorXd>(inputs[0]);

    // compute the quadratic function
    outputs[0] = (in.transpose()*Q*in + a.transpose()*in + b) (0);
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
  EXPECT_EQ(poly->numInputs, 1);
  EXPECT_EQ(poly->numOutputs, 1);

  // input
  const Eigen::VectorXd in = Eigen::VectorXd::Random(N);

  auto result = poly->Evaluate(in);

  EXPECT_DOUBLE_EQ(boost::any_cast<double>(result[0]), (in.transpose()*Q*in + a.transpose()*in + b) (0));
}
