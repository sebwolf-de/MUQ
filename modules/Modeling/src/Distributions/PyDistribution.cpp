#include "MUQ/Modeling/Distributions/PyDistribution.h"

using namespace muq::Modeling;

PyDistribution::PyDistribution(unsigned int varSizeIn, Eigen::VectorXi const& hyperSizesIn) : Distribution(varSizeIn, hyperSizesIn) {}

Eigen::VectorXd PyDistribution::SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  return SampleImpl(ToStdVec(inputs));
}

double PyDistribution::LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  return LogDensityImpl(ToStdVec(inputs));
}

std::vector<Eigen::VectorXd> PyDistribution::ToStdVec(ref_vector<Eigen::VectorXd> const& input) {
  std::vector<Eigen::VectorXd> newIns(input.size());

  for (int i=0; i<input.size(); ++i)
    newIns.at(i) = input.at(i).get();

  return newIns;
}


Eigen::MatrixXd PyGaussianBase::ApplyCovariance(Eigen::Ref<const Eigen::MatrixXd> const& x) const
{
  return ApplyCovariance(Eigen::MatrixXd(x));
};

Eigen::MatrixXd PyGaussianBase::ApplyPrecision(Eigen::Ref<const Eigen::MatrixXd> const& x) const
{
  return ApplyPrecision(Eigen::MatrixXd(x));
}

Eigen::MatrixXd PyGaussianBase::ApplyCovSqrt(Eigen::Ref<const Eigen::MatrixXd> const& x) const
{
  return ApplyCovSqrt(Eigen::MatrixXd(x));
}

Eigen::MatrixXd PyGaussianBase::ApplyPrecSqrt(Eigen::Ref<const Eigen::MatrixXd> const& x) const
{
  return ApplyPrecSqrt(Eigen::MatrixXd(x));
}

void PyGaussianBase::ResetHyperparameters(ref_vector<Eigen::VectorXd> const& params)
{
  ResetHyperparameters(PyDistribution::ToStdVec(params));
}
