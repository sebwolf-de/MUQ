#include "MUQ/Modeling/Distributions/Density.h"

using namespace muq::Modeling;


Density::Density(std::shared_ptr<Distribution> distIn) : Distribution(distIn->varSize, distIn->hyperSizes),
                                                         ModPiece(GetInputSizes(distIn),
                                                                  Eigen::VectorXi::Ones(1)),
                                                         dist(distIn)
{
  assert(dist);
}

void Density::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
  outputs.resize(1);
  outputs.at(0) = dist->LogDensity(inputs)*Eigen::VectorXd(1);
}


void Density::GradientImpl(unsigned int                const  outputDimWrt,
                           unsigned int                const  inputDimWrt,
                           ref_vector<Eigen::VectorXd> const& input,
                           Eigen::VectorXd             const& sensitivity)
{
  gradient = sensitivity(0)*dist->GradLogDensity(inputDimWrt, input);
}

void Density::JacobianImpl(unsigned int                const  outputDimWrt,
                           unsigned int                const  inputDimWrt,
                           ref_vector<Eigen::VectorXd> const& input)
{
  jacobian = dist->GradLogDensity(inputDimWrt, input).transpose();
}

void Density::ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                 unsigned int                const  inputDimWrt,
                                 ref_vector<Eigen::VectorXd> const& input,
                                 Eigen::VectorXd             const& vec)
{
  jacobianAction = dist->GradLogDensity(inputDimWrt, input).transpose() * vec;
}

double Density::LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
  return dist->LogDensity(inputs);
};

Eigen::VectorXd Density::GradLogDensityImpl(unsigned int wrt, ref_vector<Eigen::VectorXd> const& inputs)
{
  return dist->GradLogDensity(wrt, inputs);
}

Eigen::VectorXd Density::SampleImpl(ref_vector<Eigen::VectorXd> const& inputs)
{
  return dist->Sample(inputs);
}

Eigen::VectorXi Density::GetInputSizes(std::shared_ptr<Distribution> distIn)
{
  Eigen::VectorXi sizes(distIn->hyperSizes.size()+1);
  sizes(0) = distIn->varSize;
  sizes.tail(distIn->hyperSizes.size()) = distIn->hyperSizes;
  return sizes;
}
