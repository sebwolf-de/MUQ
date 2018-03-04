#ifndef KERNELIMPL_H
#define KERNELIMPL_H

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"

#include <stan/math/fwd/scal.hpp>


namespace muq
{
namespace Approximation
{

/**

\class KernelImpl
\ingroup CovarianceKernels
\brief Base class in CRTP pattern for covariance kernels
\details This class provides common functionality (such as computing Covariance matrices) for all covariance kernels.  It uses the curiously recurring template pattern and requires that child classes implement the following functions
- void EvaluateImpl(VectorType1, VectorType2, MatType)
- void GetDerivative(VectorType1, VectorType2, MatType)
- Eigen::VectorXd GetParams()
- void SetParams(Eigen::VectorXd)
*/
template<typename ChildType>
class KernelImpl : public KernelBase
{


public:

    KernelImpl(unsigned inputDimIn,
               unsigned coDimIn,
               unsigned numParamsIn) : KernelBase(inputDimIn, coDimIn, numParamsIn){};

    KernelImpl(unsigned              inputDimIn,
               std::vector<unsigned> dimIndsIn,
               unsigned              coDimIn,
               unsigned              numParamsIn) : KernelBase(inputDimIn, dimIndsIn, coDimIn, numParamsIn){};


    virtual ~KernelImpl(){};

    virtual std::shared_ptr<KernelBase> Clone() const override
    {
      return std::make_shared<ChildType>(static_cast<ChildType const &>(*this));
    }

    virtual void FillBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                           Eigen::Ref<const Eigen::VectorXd> const& x2,
                           Eigen::Ref<const Eigen::VectorXd> const& params,
                           Eigen::Ref<Eigen::MatrixXd>              block) const override
    {
      static_cast<const ChildType*>(this)->FillBlockImpl(x1,x2,params,block);
    };


    // virtual void FillPosDerivBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
    //                                Eigen::Ref<const Eigen::VectorXd> const& x2,
    //                                Eigen::Ref<const Eigen::VectorXd> const& params,
    //                                std::vector<unsigned>             const& wrts,
    //                                Eigen::Ref<Eigen::MatrixXd>              block) const override
    // {
    //   FillPosDerivBlockImpl(x1,x2,params,wrts,block);
    // }
    //
    // template<typename ScalarType>
    // virtual void FillPosDerivBlockImpl(Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>  const& x1,
    //                                    Eigen::Ref<const Eigen::VectorXd>             const& x2,
    //                                    Eigen::Ref<const Eigen::VectorXd> const&             params,
    //                                    std::vector<int>                              const& wrts,
    //                                    Eigen::Ref<Eigen::MatrixXd>                          block) const
    // {
    //   if(dimWrts.size()==0){
    //
    //     ChildType::FillBlockImpl(x1, x2, params, block);
    //
    //   }else{
    //
    //     Eigen::Matrix<stan::math::fvar<ScalarType>, Eigen::Dynamic, 1> x1_temp(x1.size());
    //     for(int i=0; i<x1.size(); ++i){
    //       x1_temp(i).val_ = x1(i);
    //       x1_temp(i).d_ = 0.0;
    //     }
    //     x1_temp(dimWrts.at(0)).d_ = 1.0;
    //
    //     std::vector<int> newWrts(dimWrts.begin()+1, dimWrts.end());
    //
    //     Eigen::Matrix<stan::math::fvar<ScalarType>, Eigen::Dynamic, Eigen::Dynamic> cov;
    //     FillPosDerivBlock(x1, x2, params, newWrts, cov);
    //
    //     for(int j=0; j<coDim.cols(); ++j){
    //       for(int i=0; i<coDim.rows(); ++i){
    //         block(i,j) = cov(i,j).d_;
    //       }
    //     }
    //
    //   }
    // };

};


}
}

#endif
