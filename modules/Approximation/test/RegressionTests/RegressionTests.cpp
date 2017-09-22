#include <gtest/gtest.h>

#include "MUQ/Approximation/Regression/Regression.h"

using namespace muq::Approximation;

class RegressionTest : public::testing::Test {
public:
  inline RegressionTest() : tol(1e-4), redtol(1e-1) {
    unsigned int Npts = 100;

    // create the linear portion
    L = Eigen::MatrixXd::Zero(2, 2);
    L << 0.5, 0.1, 0.2, 0.8;

    /*// resize the input and output points
    inputPts.resize(2,100);
    inputPts << 0.958868,0.806733,0.333761,-0.672064,0.777898,0.299414,0.258959,0.40124,-0.342446,-0.851678,-0.552687,0.0213719,-0.439916,0.438537,-0.0570331,0.888636,-0.327298,-0.130973,-0.310114,0.666487,0.350952,-0.0361284,0.424175,0.243646,-0.172033,0.347873,-0.305768,0.218212,0.461459,0.480877,0.841829,0.306261,0.0648819,-0.479006 ,0.37225,-0.777449,0.153381,0.333113,0.551534,-0.340716,0.968726,0.654782,-0.623598,0.917274,0.529743,-0.757714,-0.232336,0.886103,0.723834,0.587314,-0.405423,0.819286,-0.00371215,-0.674487,0.729158,-0.0726757,-0.00804521,-0.639158,0.455101,0.206218,0.676267,-0.643585,-0.00294904,-0.723523,-0.350386,0.816969,0.673656,-0.00785116,-0.211346,0.217766,-0.69754,-0.784303,-0.272803,-0.337228,-0.145345,0.167141,0.317493,-0.0251463,0.766073,0.0354294,0.115121,0.659878,-0.511346,0.45872,0.969689,0.795121,-0.178424,0.566564,-0.412644,0.73107,-0.901675,0.972934,-0.578234,0.730362,-0.800881,-0.396474,0.618191,-0.896983,-0.0845688,0.384153,
    0.487622,0.967191,-0.00548296,0.660024,-0.846011,-0.503912,-0.541726,-0.366266,-0.537144,0.266144,0.302264,0.942931,0.0922138,-0.773439 ,0.18508,-0.0981649,0.695369,-0.993537,0.196963,-0.532217,-0.0340994,-0.390089,-0.634888,-0.918271,0.391968 ,0.27528,-0.630755,0.254316,-0.343251,-0.595574,0.369513,-0.485469,-0.824713,0.754768,-0.81252,-0.276798,0.186423,-0.422444,-0.423241,-0.620498,-0.992843,-0.337042,-0.127006,0.837861,0.398151,0.371572,0.548547,0.832546,-0.592904,0.0960841,0.809865,0.747958,0.152399,-0.452178,-0.0152024,0.697884,-0.417893,0.368357,-0.721884,-0.0151566,0.448504,-0.556069,-0.757482,-0.279115,0.863791,0.244191,0.636255,-0.330057,0.317662,-0.482188,-0.85491,0.294415,-0.423461,-0.817703,0.868989,-0.469077,0.523556,-0.685456,0.251331,-0.584313,-0.147601,-0.211223,-0.347973,0.277308,-0.323514,-0.727851,-0.989183,0.548772,-0.770664,0.442012,-0.10179,0.415818,-0.052212,-0.812161,-0.234208,0.31424,-0.736596,-0.893155,0.561737,-0.11488;*/

    // generate the input points
    ins.resize(100, Eigen::Vector2d::Constant(std::numeric_limits<double>::quiet_NaN()));
    for( auto it=ins.begin(); it!=ins.end(); ++it ) { *it = Eigen::Vector2d::Random(); }

    // generate the output points
    outs.resize(ins.size(), Eigen::Vector2d::Constant(std::numeric_limits<double>::quiet_NaN()));
    for( unsigned int i=0; i<ins.size(); ++i ) {
      const Eigen::VectorXd temp = L*ins[i];
      outs[i](0) = std::exp(temp(0));
      outs[i](1) = std::cos(temp(1));
    }

    /*// generate the points
    outputPts = Eigen::MatrixXd::Zero(2, Npts);
    for (unsigned int i = 0; i < Npts; ++i) {
      Eigen::VectorXd temp = L * inputPts.col(i);
      outputPts(0, i) = exp(temp(0));
      outputPts(1, i) = cos(temp(1));
      }*/

    // set the test input location
    testIn  = Eigen::VectorXd::Random(2);
    trueOut = Eigen::VectorXd::Zero(2);
    trueJac = Eigen::MatrixXd::Zero(2, 2);

    // compute the true output
    Eigen::VectorXd temp = L * testIn;
    trueOut(0) = exp(temp(0));
    trueOut(1) = cos(temp(1));

    // compute the true jacobian
    trueJac(0, 0) = trueOut(0) * L(0, 0);
    trueJac(0, 1) = trueOut(0) * L(0, 1);
    trueJac(1, 0) = -1.0 * L(1, 0) * sin(temp(1));
    trueJac(1, 1) = -1.0 * L(1, 1) * sin(temp(1));
  }

  inline virtual ~RegressionTest() {}

  /// A matrix holding the input points.
  std::vector<Eigen::Vector2d> ins;
  //Eigen::MatrixXd inputPts;

  /// A matrix holding the output points.
  std::vector<Eigen::Vector2d> outs;
  //Eigen::MatrixXd outputPts;

  /// The linear part of the true forward model.
  Eigen::MatrixXd L;

  /// The locaiton to evaluate the approximation for testing.
  Eigen::VectorXd testIn;

  /// Store the true output at the test input. 
  Eigen::VectorXd trueOut;

  /// Store the true jacobian at the test input
  Eigen::MatrixXd trueJac;

  /// Define the tolerance allowed in the approximation at testIn.
  const double tol;

  /// Define a tolerance for reduced order regressions.
  const double redtol;
};

TEST_F(RegressionTest, Fit) {
  // create the regression 
  auto reg = std::make_shared<Regression>();

  // fit the polynomial coefficients
  reg->Fit<Eigen::Vector2d>(ins, outs);
  std::cout << "!!!!!!!!!!!!!!!!!!!" << std::endl;
  reg->Fit<Eigen::Vector2d>(ins, outs, Eigen::Vector2d::Ones());
}

