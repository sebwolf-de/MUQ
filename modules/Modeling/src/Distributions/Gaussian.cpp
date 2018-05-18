#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

Gaussian::Gaussian(unsigned int dim,
                   InputMask    extraInputs) : Distribution(dim, GetExtraSizes(dim, extraInputs)),
                                               mode(ModeFromExtras(extraInputs)),
                                               inputTypes(extraInputs),
                                               mean(Eigen::VectorXd::Zero(dim)),
                                               covPrec(Eigen::VectorXd::Ones(dim))
{
  ComputeNormalization();
}

Gaussian::Gaussian(Eigen::VectorXd const& muIn,
                   InputMask              extraInputs) : Distribution(muIn.size(), GetExtraSizes(muIn.size(), extraInputs)),
                                                         mode(ModeFromExtras(extraInputs)),
                                                         inputTypes(extraInputs),
                                                         mean(muIn),
                                                         covPrec(Eigen::VectorXd::Ones(muIn.size()))
{
  ComputeNormalization();
}


Gaussian::Gaussian(Eigen::VectorXd const& muIn,
                   Eigen::MatrixXd const& obj,
                   Gaussian::Mode         modeIn,
                   InputMask              extraInputs) : Distribution(muIn.size(), GetExtraSizes(muIn.size(), extraInputs)),
                                                         mode(modeIn),
                                                         inputTypes(extraInputs),
                                                         mean(muIn),
                                                         covPrec(obj)
{
  CheckInputTypes(extraInputs, mode);
  assert(mean.rows()==covPrec.rows());
  if(covPrec.cols()>1)
    assert(mean.rows()==covPrec.cols());

  if(covPrec.cols()>1)
    sqrtCovPrec = covPrec.selfadjointView<Eigen::Lower>().llt();

  ComputeNormalization();
}

void Gaussian::CheckInputTypes(InputMask extraInputs, Mode mode)
{
  if( (((extraInputs & ExtraInputs::DiagCovariance)>0) || ((extraInputs & ExtraInputs::FullCovariance)>0)) && (mode == Mode::Precision))
    throw std::logic_error("Extra arguments passed to Gaussian constructor do not match the covariance mode.");
  if( (((extraInputs & ExtraInputs::DiagPrecision)>0) || ((extraInputs & ExtraInputs::FullPrecision)>0)) && (mode == Mode::Covariance))
    throw std::logic_error("Extra arguments passed to Gaussian constructor do not match the covariance mode.");
}

Gaussian::Mode Gaussian::ModeFromExtras(InputMask extraInputs)
{
  if( ((extraInputs & ExtraInputs::DiagCovariance)>0) || ((extraInputs & ExtraInputs::FullCovariance)>0)){
    return Gaussian::Mode::Covariance;
  }else{
    return Gaussian::Mode::Precision;
  }
}

Eigen::VectorXi Gaussian::GetExtraSizes(unsigned dim, InputMask extraInputs)
{
  Eigen::VectorXi output(2);
  int numExtras = 0;

  if((extraInputs & ExtraInputs::Mean)>0){
    output(0) = dim;
    numExtras++;
  }

  if( ((extraInputs & ExtraInputs::DiagCovariance)>0) || ((extraInputs & ExtraInputs::DiagPrecision)>0)){
    output(numExtras) = dim;
    numExtras++;
  }


  if( ((extraInputs & ExtraInputs::FullCovariance)>0) || ((extraInputs & ExtraInputs::FullPrecision)>0)){
    assert(numExtras<2);
    output(numExtras) = dim*dim;
    numExtras++;
  }

  return output.head(numExtras);
}


Eigen::MatrixXd Gaussian::GetCovariance() const
{
  if(mode==Mode::Covariance){
    if(covPrec.cols()==1){
      return covPrec.col(0).asDiagonal();
    }else{
      return covPrec;
    }
  }else{
    if(covPrec.cols()==1){
      return covPrec.col(0).array().inverse().matrix().asDiagonal();
    }else{
      return covPrec.selfadjointView<Eigen::Lower>().llt().solve(Eigen::MatrixXd::Identity(mean.rows(),mean.rows()));
    }
  }
}

Eigen::MatrixXd Gaussian::GetPrecision() const
{
  if(mode==Mode::Precision){
    if(covPrec.cols()==1){
      return covPrec.col(0).asDiagonal();
    }else{
      return covPrec;
    }
  }else{
    if(covPrec.cols()==1){
      return covPrec.col(0).array().inverse().matrix().asDiagonal();
    }else{
      return covPrec.selfadjointView<Eigen::Lower>().llt().solve(Eigen::MatrixXd::Identity(mean.rows(),mean.rows()));
    }
  }
}


void Gaussian::ComputeNormalization() {
  logNormalization = -0.5*Dimension()*std::log(2.0*M_PI);

  if( mode==Gaussian::Mode::Covariance ) {
    if(covPrec.cols()==1){
      logNormalization -= 0.5 * covPrec.array().log().sum();
    }else{
      logNormalization -= std::log(sqrtCovPrec.matrixL().determinant());
    }
  } else if( mode==Gaussian::Mode::Precision ) {
    if(covPrec.cols()==1){
      logNormalization += 0.5 * covPrec.array().log().sum();
    }else{
      logNormalization += std::log(sqrtCovPrec.matrixL().determinant());
    }
  }
}

void Gaussian::ResetHyperparameters(ref_vector<Eigen::VectorXd> const& inputs)
{
  unsigned currInd = 0;

  if((inputTypes & ExtraInputs::Mean)>0){
    assert(inputs.at(currInd).get().size() == mean.size());
    mean = inputs.at(currInd);
    currInd++;
  }

  if(((inputTypes & ExtraInputs::DiagCovariance)>0) || ((inputTypes & ExtraInputs::DiagPrecision)>0)){
    assert(inputs.at(currInd).get().size() == covPrec.rows());
    covPrec = inputs.at(currInd).get();
  }else if(((inputTypes & ExtraInputs::FullCovariance)>0) || ((inputTypes & ExtraInputs::FullPrecision)>0)){
    assert(inputs.at(currInd).get().size() == covPrec.rows()*covPrec.cols());
    Eigen::Map<const Eigen::MatrixXd> mat(inputs.at(currInd).get().data(), covPrec.rows(), covPrec.cols());
    covPrec = mat;
    sqrtCovPrec = covPrec.selfadjointView<Eigen::Lower>().llt();
  }

  ComputeNormalization();
}

double Gaussian::LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) {

  ResetHyperparameters(ref_vector<Eigen::VectorXd>(inputs.begin()+1, inputs.end()));

  Eigen::VectorXd delta = inputs.at(0).get() - mean;

  switch( mode ) {
    case Gaussian::Mode::Covariance: {
      if(covPrec.cols()==1){
        return logNormalization - 0.5 * delta.dot( (covPrec.col(0).array().inverse() * delta.array()).matrix());
      }else{
        return logNormalization - 0.5 * delta.dot( sqrtCovPrec.solve(delta));
      }
    }
    case Gaussian::Mode::Precision: {
      if(covPrec.cols()==1){
        return logNormalization - 0.5 * delta.dot( (covPrec.col(0).array() * delta.array()).matrix());
      }else{
        return logNormalization - 0.5 * delta.dot( covPrec.selfadjointView<Eigen::Lower>() * delta );
      }
    }
    default: {
      assert(false);
      return -std::numeric_limits<double>::infinity();
    }
  }
}

Eigen::VectorXd Gaussian::SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) {

  ResetHyperparameters(ref_vector<Eigen::VectorXd>(inputs.begin(), inputs.end()));

  Eigen::VectorXd z = RandomGenerator::GetNormal(mean.rows());

  switch( mode ) {
  case Gaussian::Mode::Covariance: {
    if(covPrec.cols()==1){
      return mean + covPrec.col(0).array().sqrt().matrix().asDiagonal() * z;
    }else{
      return mean + sqrtCovPrec.matrixL() * z;
    }
  }
  case Gaussian::Mode::Precision: {
    if(covPrec.cols()==1){
      return mean + covPrec.col(0).array().inverse().sqrt().matrix().asDiagonal() * z;
    }else{
      return mean + sqrtCovPrec.matrixL().solve( z );
    }
  }
  default: {
    assert(false);
    return Eigen::VectorXd();
  }
  }
}

unsigned int Gaussian::Dimension() const {
  return mean.rows();
}

void Gaussian::SetMean(Eigen::VectorXd const& newMu)
{
  assert(newMu.rows() == mean.rows());

  mean = newMu;
}

void Gaussian::SetCovariance(Eigen::MatrixXd const& newCov) {

  mode = Gaussian::Mode::Covariance;

  assert(newCov.rows() == mean.rows());

  if(newCov.cols()>1)
    assert(newCov.cols() == mean.rows());

  covPrec = newCov;
  sqrtCovPrec = covPrec.selfadjointView<Eigen::Lower>().llt();

  // recompute the scaling constant
  ComputeNormalization();
}

void Gaussian::SetPrecision(Eigen::MatrixXd const& newPrec) {

  mode = Gaussian::Mode::Precision;

  assert(newPrec.rows() == mean.rows());

  if(newPrec.cols()>1)
    assert(newPrec.cols() == mean.rows());

  covPrec = newPrec;
  sqrtCovPrec = covPrec.selfadjointView<Eigen::Lower>().llt();

  // recompute the scaling constant
  ComputeNormalization();
}
