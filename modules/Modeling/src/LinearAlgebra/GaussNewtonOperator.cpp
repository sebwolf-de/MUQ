#include "MUQ/Modeling/LinearAlgebra/GaussNewtonOperator.h"

#include "MUQ/Modeling/Distributions/Density.h"

using namespace muq::Modeling;


GaussNewtonOperator::GaussNewtonOperator(std::shared_ptr<ModPiece>     const& forwardModelIn,
                                         std::shared_ptr<GaussianBase> const& noiseModelIn,
                                         std::vector<Eigen::VectorXd>  const& inputsIn,
                                         unsigned int                         inWrtIn) : LinearOperator(forwardModelIn->inputSizes(inWrt), forwardModelIn->inputSizes(inWrt)),
                                                                                         forwardModel(forwardModelIn),
                                                                                         noiseModel(noiseModelIn),
                                                                                         inputs(inputsIn),
                                                                                         inWrt(inWrtIn)
{
  assert(forwardModel->outputSizes.size()==1);
  assert(forwardModel->outputSizes(0)==noiseModelIn->AsDensity()->inputSizes(0));
  assert(forwardModel->inputSizes.size()>inWrtIn);
}


Eigen::MatrixXd GaussNewtonOperator::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
  Eigen::MatrixXd output(rows(),x.cols());
  for(unsigned int i=0; i<x.cols(); ++i){
    Eigen::VectorXd temp =forwardModel->ApplyJacobian(0,inWrt,inputs,x.col(i).eval());
    temp = noiseModel->ApplyPrecision(temp);
    output.col(i) = forwardModel->Gradient(0,inWrt,inputs,temp);
  }
  return output;
}

Eigen::MatrixXd GaussNewtonOperator::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
  return Apply(x);
}
