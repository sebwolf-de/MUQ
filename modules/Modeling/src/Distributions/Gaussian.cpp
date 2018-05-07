#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

Gaussian::Gaussian(boost::any               const& obj,
                   Gaussian::Mode           const  mode,
                   Gaussian::ExpectedInputs const  expected) :
      Distribution(),
      mode(mode),
      dim(algebra->Size(obj,0)),
      mean(mode==Gaussian::Mode::Mean? obj : boost::none),
      cov(mode==Gaussian::Mode::Covariance? SaveCovPrec(obj) : boost::none),
      prec(mode==Gaussian::Mode::Precision? SaveCovPrec(obj) : boost::none)
{
  SetInputSize(0,dim);
  SetInputTypes(expected);
  ComputeScalingConstant();
}

Gaussian::Gaussian(boost::any               const& mean,
                   boost::any               const& obj,
                   Gaussian::Mode           const  mode,
                   Gaussian::ExpectedInputs const  expected) :
      Distribution(),
      mode(mode),
      dim(algebra->Size(obj,0)),
      mean(mean),
      cov(mode==Gaussian::Mode::Covariance? SaveCovPrec(obj) : boost::none),
      prec(mode==Gaussian::Mode::Precision? SaveCovPrec(obj) : boost::none)
{
  SetInputTypes(expected);
  ComputeScalingConstant();
}

Gaussian::~Gaussian() {}


void Gaussian::SetInputTypes(Gaussian::ExpectedInputs const  expected)
{
  if(expected == ExpectedInputs::LocationOnly){
    numInputs = 2;

  }else if(expected == ExpectedInputs::LocationMean){
    numInputs = 3;

  }else if(expected == ExpectedInputs::All){
    numInputs = 4;

  }else if(expected == ExpectedInputs::Flexible){
    numInputs = -1;
  }
}


void Gaussian::ComputeScalingConstant() {
  scalingConstant = dim*std::log(2.0*M_PI);

  if( mode==Gaussian::Mode::Covariance ) {
    scalingConstant += algebra->LogDeterminate(*cov);
  } else if( mode==Gaussian::Mode::Precision ) {
    scalingConstant -= algebra->LogDeterminate(*prec);
  }
}

boost::any Gaussian::SaveCovPrec(boost::any const& in) {
  // for Eigen matrices, we can pre compute the cholesky
  if( in.type()==typeid(Eigen::MatrixXd) ) {
    Eigen::LLT<Eigen::MatrixXd> chol;
    chol.compute(boost::any_cast<Eigen::MatrixXd const&>(in));

    return chol;
  }

  return in;
}

void Gaussian::ResetHyperparameters(ref_vector<boost::any> const& hyperparas) {
  // reset the hyper parameters
  for( auto hp : hyperparas ) {
    // get the hyperparameter
    const std::pair<boost::any, Gaussian::Mode>& hyperpara = boost::any_cast<std::pair<boost::any, Gaussian::Mode> >(hp);

    switch( hyperpara.second ) {
    case Gaussian::Mode::Mean: {
      // reset the mean
      mean = hyperpara.first;
      break;
    }
    case Gaussian::Mode::Covariance: {
      SetCovariance(hyperpara.first);
      break;
    }
    case Gaussian::Mode::Precision: {
      SetPrecision(hyperpara.first);
      break;
    }
    default: {
      // something is wrong
      assert(false);
      break;
    }
    }
  }
}

double Gaussian::LogDensityImpl(ref_vector<boost::any> const& inputs) {
  ResetHyperparameters(ref_vector<boost::any>(inputs.begin()+1, inputs.end()));

  boost::any delta = inputs[0].get();
  if( mean ) { delta = algebra->Subtract(*mean, delta); }

  switch( mode ) {
  case Gaussian::Mode::Covariance: {
    return -0.5*(scalingConstant+algebra->InnerProduct(delta, algebra->ApplyInverse(*cov, delta)));
  }
  case Gaussian::Mode::Precision: {
    return -0.5*(scalingConstant+algebra->InnerProduct(delta, algebra->Apply(*prec, delta)));
  }
  case Gaussian::Mode::Mean: {
    return -0.5*(scalingConstant+algebra->InnerProduct(delta, delta));
  }
  default: {
    assert(false);
    return 0.0;
  }
  }
}

boost::any Gaussian::SampleImpl(ref_vector<boost::any> const& inputs) {
  ResetHyperparameters(ref_vector<boost::any>(inputs.begin(), inputs.end()));

  // make sure the dimension is nonzero
  assert(dim>0);

  boost::any stdnrm;
  if( dim==1 ) {
    stdnrm = (double)RandomGenerator::GetNormal();
  } else {
    stdnrm = (Eigen::VectorXd)RandomGenerator::GetNormal(dim);
  }

  switch( mode ) {
  case Gaussian::Mode::Covariance: {
    if( !covSqrt ) { covSqrt = algebra->SquareRoot(*cov); }
    return mean? algebra->Add(algebra->Apply(*covSqrt, stdnrm), *mean) : algebra->Apply(*covSqrt, stdnrm);
  }
  case Gaussian::Mode::Precision: {
    if( !precSqrt ) { precSqrt = algebra->SquareRoot(*prec); }

    return mean? algebra->Add(algebra->ApplyInverse(*precSqrt, stdnrm), *mean) : algebra->ApplyInverse(*precSqrt, stdnrm);
  }
  case Gaussian::Mode::Mean: {
    return algebra->Add(stdnrm, *mean);
  }
  default: {
    assert(false);
    return boost::none;
  }
  }
}

unsigned int Gaussian::Dimension() const {
  return dim;
}

boost::any Gaussian::GetCovariance() const {
  assert(cov);
  return *cov;
}

void Gaussian::SetCovariance(boost::any const& newcov) {
  // reset the mode
  mode = Gaussian::Mode::Covariance;

  // reset the covaraince and the precision
  prec = boost::none;
  precSqrt = boost::none;
  cov = SaveCovPrec(newcov);
  covSqrt = boost::none;

  // recompute the scaling constant
  ComputeScalingConstant();
}

void Gaussian::SetPrecision(boost::any const& newprec) {
  // reset the mode
  mode = Gaussian::Mode::Precision;

  // reset the covaraince and the precision
  prec = SaveCovPrec(newprec);
  precSqrt = boost::none;
  cov = boost::none;
  covSqrt = boost::none;

  // recompute the scaling constant
  ComputeScalingConstant();
}
