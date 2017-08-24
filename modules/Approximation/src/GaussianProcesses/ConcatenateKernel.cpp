#include "MUQ/Approximation/GaussianProcesses/ConcatenateKernel.h"

ConcatenateKernel::ConcatenateKernel(std::vector<std::shared_ptr<KernelBase>> const& kernelsIn) : KernelBase(kernelsIn.at(0)->inputDim,
                                                                                                             GetNumCodim(kernelsIn),
                                                                                                             GetNumParams(kernelsIn)),
                                                                                                  kernels(kernelsIn)
{};

Eigen::MatrixXd ConcatenateKernel::Evaluate(Eigen::VectorXd const& x1, Eigen::VectorXd const& x2) const override
{
  Eigen::MatrixXd output = Eigen::MatrixXd::Zero(coDim,coDim);

  int currRow = 0;
  for(int i=0; i<kernels.size(); ++i){
    int currCoDim = kernels.at(i)->coDim;
    output.block(currRow, currRow, currCoDim, currCoDim) = kernel.at(i)->Evaluate(x1,x2);
    currRow += currCoDim;
  }

  return output;
}

virtual void ConcatenateKernel::FillDerivativeMatrix(Eigen::MatrixXd             const& xs,
                                                        unsigned                           wrt,
                                                        Eigen::Ref<Eigen::MatrixXd>        derivs) const override
{
  assert(wrt < this->numParams);
  
  // Figure out which kernel we're looking at
  int wrtKernel = 0;
  int paramCum = 0;
  int coCum = 0;
  for(int i=0; i<kernels.size(); ++i)
  {
    paramCum += kernels.at(i)->numParams;
    if(wrt<paramCum){
      wrtKernel = i;
      break;
    }
    coCum += kernels.at(i)->coDim;
  }
  
  // Initialize the derivative matrix to 0
  for(int j=0; j<derivs.cols(); ++j)
  {
    for(int i=0; i<derivs.rows(); ++i)
      derivs(i,j) = 0.0;
  }

  kernels.at(wrtKernel)->FillDerivativeMatrix(derivs.block(coCum, coCum, kernels.at(wrtKernel)->coDim, kernels.at(wrtKernel)->coDim));
}


int ConcatenateKernel::GetNumCodim(std::vector<std::shared_ptr<KernelBase>> const& kernelsIn)
{
  int count = 0;
  for(auto kptr : kernelsIn)
    count += kptr->coDim;
  return count;
}

int ConcatenateKernel::GetNumParams(std::vector<std::shared_ptr<KernelBase>> const& kernelsIn)
{
  int count = 0;
  for(auto kptr : kernelsIn)
    count += kptr->numParams;
  return count;
}

