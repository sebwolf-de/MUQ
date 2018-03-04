#ifndef COREGIONALKERNEL_H
#define COREGIONALKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"


namespace muq
{
namespace Approximation
{

/**

@class CoregionalKernel
@ingroup CovarianceKernels
@brief This kernel supports coregionalization for modeling vector-valued Gaussian processes.
@details Let \f$v_i\f$ and \f$\lambda_i\f$ be the \f$i^{th}\f$ eigenvector and eigenvalue of an \f$N\times N\f$ covariance matrix \f$\Sigma\f$.  Now, let define the matrix \f$\Phi\f$ as
\f[
\Phi = \left[\sqrt{\lambda_1} v_1, \sqrt{\lambda_2} v_2, \ldots, \sqrt{\lambda_N} v_N \right]
\f]
The CoregionalKernels defined by this class take the form
\f[
k(x,x^\prime) = \Phi \left[\begin{array}{cccc} k_1(x,x^\prime) & 0 & & 0\\ 0 & k_2(x,x^\prime) & \ddots & \vdots\\ \vdots & & \ddots  & \\ 0 & \cdots & 0 & k_N(x,x^\prime) \end{array}\right] \Phi^T.
\f]
Notice that the marginal covariance (i.e., \f$k(x,x)\f$) is \f$\Sigma\f$.
 */
class CoregionalKernel : public KernelBase
{

public:

    CoregionalKernel(unsigned                                        dim,
		                 Eigen::MatrixXd                          const& Gamma,
		                 std::vector<std::shared_ptr<KernelBase>> const& kernelsIn);

    virtual ~CoregionalKernel(){};

    virtual void FillBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                           Eigen::Ref<const Eigen::VectorXd> const& x2,
                           Eigen::Ref<const Eigen::VectorXd> const& params,
                           Eigen::Ref<Eigen::MatrixXd>              block) const override
    {
	    Eigen::VectorXd rhoVec(coDim);

      unsigned currInd = 0;
      for(int i=0; i<coDim; ++i){
        kernels.at(i)->SetParams(params.segment(currInd, kernels.at(i)->numParams));
        rhoVec(i) = kernels.at(i)->Evaluate(x1,x2)(0,0);
        currInd += kernels.at(i)->numParams;
      }

	     block = A * rhoVec.asDiagonal() * A.transpose();
    }
  //
  //   template<typename VecType1, typename VecType2, typename MatrixType>
  //   inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatrixType & derivs) const
  //   {
	// unsigned cumParams = 0;
  //       unsigned kernelInd = 0;
  //       Eigen::MatrixXd kernelDerivs;
  //
  //       for(unsigned i=0; i<this->coDim; ++i)
  //       {
  //           if(wrt < cumParams + kernels.at(i)->numParams)
  //           {
  //               kernelDerivs = kernels.at(i)->GetDerivative(x1, x2, wrt-cumParams);
  //               kernelInd = i;
  //               break;
  //           }
  //           else
  //           {
  //               cumParams += kernels.at(i)->numParams;
  //           }
  //       }
  //
  //       assert(kernelDerivs.rows()==1);
  //       assert(kernelDerivs.cols()==1);
  //
	// derivs = Eigen::MatrixXd(A.col(kernelInd) * kernelDerivs(0,0) * A.col(kernelInd).transpose());
  //   }

  virtual std::shared_ptr<KernelBase> Clone() const override{return std::make_shared<CoregionalKernel>(*this);};
  
private:

    Eigen::MatrixXd A;
    std::vector<std::shared_ptr<KernelBase>> kernels;

    static unsigned GetNumParams(std::vector<std::shared_ptr<KernelBase>> const& kernelsIn)
    {
        unsigned allParams = 0;
        for(unsigned i = 0; i<kernelsIn.size(); ++i)
            allParams += kernelsIn.at(i)->numParams;
        return allParams;
    };

};


template<class KernelType1, class... KTypes>
void FillKernelVector(std::vector<std::shared_ptr<KernelBase>> & vec, const KernelType1& kernel1, const KTypes&... kernels)
{
    vec.push_back(kernel1.Clone());
    FillKernelVector(vec, kernels...);
}

template<class KernelType1>
void FillKernelVector(std::vector<std::shared_ptr<KernelBase>> & vec, const KernelType1& kernel1)
{
    KernelBase const& baseRef = kernel1;
    vec.push_back(baseRef.Clone());
}


/**

Create a coregional kernel
@param[in] coCov The joint covariance of the vector-valued process
@param[in] kernels Covariance kernels for the principal components of the covariance.

*/
template<class KernelType1, class KernelType2, typename = typename std::enable_if<std::is_base_of<KernelBase, KernelType2>::value>, class... KTypes>
CoregionalKernel CoregionTie(Eigen::MatrixXd const& coCov, const KernelType1& kernel1, const KernelType2& kernel2, const KTypes&... kernels)
{
    std::vector<std::shared_ptr<KernelBase>> vec;
    FillKernelVector(vec, kernel1, kernel2, kernels...);

    return CoregionalKernel(kernel1.inputDim, coCov, vec);

}
template<class KernelType1>
CoregionalKernel CoregionTie(Eigen::MatrixXd const& coCov, const KernelType1& kernel1, int numRepeat)
{
    std::vector<std::shared_ptr<KernelBase>> vec(numRepeat);
    KernelBase const& baseRef = kernel1;
    for(int i=0; i<numRepeat; ++i)
        vec.at(i) = baseRef.Clone();

    return CoregionalKernel(kernel1.inputDim, coCov, vec);
}

}
}


#endif
