#ifndef PRODUCTKERNEL_H
#define PRODUCTKERNEL_H

#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"
#include "MUQ/Approximation/GaussianProcesses/PeriodicKernel.h"
#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"

#include "MUQ/Modeling/LinearSDE.h"

#include "MUQ/Utilities/LinearAlgebra/KroneckerProductOperator.h"
#include "MUQ/Utilities/LinearAlgebra/BlockDiagonalOperator.h"
#include "MUQ/Utilities/LinearAlgebra/BlockRowOperator.h"
#include "MUQ/Utilities/LinearAlgebra/SumOperator.h"
#include "MUQ/Utilities/LinearAlgebra/IdentityOperator.h"



#include "MUQ/Utilities/Exceptions.h"

namespace muq
{
namespace Approximation
{

/**

@class ProductKernel
#ingroup CovarianceKernels
\f[
k(x,y) = k_1(x,y)*k_2(x,y)
\f]

 */
class ProductKernel : public KernelBase
{

public:
    ProductKernel(std::shared_ptr<KernelBase> kernel1In,
                  std::shared_ptr<KernelBase> kernel2In);

    virtual ~ProductKernel(){};

    virtual void FillBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                           Eigen::Ref<const Eigen::VectorXd> const& x2,
                           Eigen::Ref<const Eigen::VectorXd> const& params,
                           Eigen::Ref<Eigen::MatrixXd>              block) const override;

  //   template<typename VecType1, typename VecType2, typename MatType>
  //   inline void GetDerivative(VecType1 const& x1, VecType2 const& x2, int wrt, MatType & derivs) const
  //   {
	// assert(wrt < this->numParams);
  //
	// Eigen::MatrixXd temp1(kernel1.coDim, kernel1.coDim);
	// Eigen::MatrixXd temp2(kernel2.coDim, kernel2.coDim);
  //
	// if(wrt < kernel1.numParams )
	// {
	//     kernel1.GetDerivative(x1, x2, wrt, temp1);
	//     kernel2.EvaluateImpl(x1,x2, temp2);
	// }
	// else
	// {
	//     kernel1.EvaluateImpl(x1,x2, temp1);
	//     kernel2.GetDerivative(x1,x2, wrt-kernel1.numParams, temp2);
	// }
  //
	// if(kernel1.coDim==kernel2.coDim)
	// {
	//     derivs = Eigen::MatrixXd(temp1.array() * temp2.array());
	// }
	// else if(kernel1.coDim==1)
	// {
	//     derivs = temp1(0,0) * temp2;
	// }
	// else if(kernel2.coDim==1)
	// {
  //           derivs = temp2(0,0) * temp1;
	// }
	// else
	// {
	//     std::cerr << "\nERROR: Something unexpected happened with the dimensions of the kernels in this product.\n";
	//     assert(false);
	// }

    //}

    virtual std::vector<std::shared_ptr<KernelBase>> GetSeperableComponents() override;

    virtual std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Utilities::LinearOperator>, Eigen::MatrixXd> GetStateSpace(boost::property_tree::ptree sdeOptions=boost::property_tree::ptree()) const override;

protected:
    std::shared_ptr<KernelBase> kernel1;
    std::shared_ptr<KernelBase> kernel2;

    std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Utilities::LinearOperator>, Eigen::MatrixXd> GetProductStateSpace(std::shared_ptr<PeriodicKernel> const& kernel1,
                                                                                                                                                 	              std::shared_ptr<KernelBase>     const& kernel2,
                                                                                                                                                                boost::property_tree::ptree sdeOptions) const;

};


}
}


#endif
