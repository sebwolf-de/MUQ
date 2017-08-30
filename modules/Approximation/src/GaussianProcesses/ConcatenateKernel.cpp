#include "MUQ/Approximation/GaussianProcesses/ConcatenateKernel.h"

using namespace muq::Approximation;


ConcatenateKernel::ConcatenateKernel(std::vector<std::shared_ptr<KernelBase>> const& kernelsIn) : KernelBase(kernelsIn.at(0)->inputDim,
                                                                                                             GetNumCodim(kernelsIn),
                                                                                                             GetNumParams(kernelsIn)),
                                                                                                  kernels(kernelsIn)
{};

Eigen::MatrixXd ConcatenateKernel::Evaluate(Eigen::Ref<const Eigen::VectorXd> const& x1, 
                                             Eigen::Ref<const Eigen::VectorXd> const& x2) const
{
  Eigen::MatrixXd output = Eigen::MatrixXd::Zero(coDim,coDim);

  int currRow = 0;
  for(int i=0; i<kernels.size(); ++i){
    int currCoDim = kernels.at(i)->coDim;
    output.block(currRow, currRow, currCoDim, currCoDim) = kernels.at(i)->Evaluate(x1,x2);
    currRow += currCoDim;
  }

  return output;
}

void ConcatenateKernel::FillCovariance(Eigen::MatrixXd             const& xs,
                                       Eigen::MatrixXd             const& ys,
                                       Eigen::Ref<Eigen::MatrixXd>        cov) const
{
  
  for(int j=0; j<cov.cols(); ++j){
    for(int i=0; i<cov.rows(); ++i)
      cov(i,j) = 0.0;
  }

  Eigen::MatrixXd tempCov;

  // Loop over all the pairs of points
  for(int xind=0; xind<xs.cols(); ++xind){
    for(int yind=0; yind<ys.cols(); ++yind){
      tempCov = Evaluate(xs.col(xind), ys.col(yind));
      
      for(int j=0; j<coDim; ++j){
        for(int i=0; i<coDim; ++i)
          cov(xind*coDim + i, yind*coDim + j) = tempCov(i,j);
      }
    }
  }
}

Eigen::MatrixXd ConcatenateKernel::GetDerivative(Eigen::Ref<const Eigen::VectorXd> const& x1, 
                                                 Eigen::Ref<const Eigen::VectorXd> const& x2, 
                                                 int                                wrt) const
{
    assert(wrt < this->numParams);
  
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(coDim, coDim);

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

    output.block(coCum, coCum, kernels.at(wrtKernel)->coDim, kernels.at(wrtKernel)->coDim) = kernels.at(wrtKernel)->GetDerivative(x1, x2, wrt + kernels.at(wrtKernel)->numParams - paramCum);

    return output;
}

void ConcatenateKernel::FillCovariance(Eigen::MatrixXd             const& xs,
                                       Eigen::Ref<Eigen::MatrixXd>        cov) const
{
  for(int j=0; j<cov.cols(); ++j){
    for(int i=0; i<cov.rows(); ++i)
      cov(i,j) = 0.0;
  }

  // Loop over all the pairs of points
  for(int xind=0; xind<xs.cols(); ++xind){
    for(int yind=0; yind<=xind; ++yind){
      cov.block(xind*coDim, yind*coDim, coDim, coDim) = Evaluate(xs.col(xind), xs.col(yind));
    }
  }
  cov.triangularView<Eigen::Upper>() = cov.triangularView<Eigen::Lower>().transpose();
}

void ConcatenateKernel::FillDerivativeMatrix(Eigen::MatrixXd             const& xs,
                                             unsigned                           wrt,
                                             Eigen::Ref<Eigen::MatrixXd>        derivs) const
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

  const int wrtCoDim = kernels.at(wrtKernel)->coDim;

  for(int j=0; j<xs.cols(); ++j){
    for(int i=0; i<xs.cols(); ++i){
      derivs.block(i*coDim+coCum, j*coDim+coCum, wrtCoDim, wrtCoDim) = kernels.at(wrtKernel)->GetDerivative(xs.col(i), xs.col(j), wrt + kernels.at(wrtKernel)->numParams - paramCum);
    
    }
  }  
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

std::shared_ptr<KernelBase> ConcatenateKernel::Clone() const
{
  std::vector<std::shared_ptr<KernelBase>> newKernels(kernels.size());
  for(int i=0; i<kernels.size(); ++i)
    newKernels.at(i) = kernels.at(i)->Clone();
 
  return std::make_shared<ConcatenateKernel>(newKernels);
}
