#ifndef DENSELINEAROPERATOR_H
#define DENSELINEAROPERATOR_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

namespace muq
{
namespace Modeling
{


/** @class DenseLinearOperator
 *  @ingroup Utilities
 *  @brief Wraps a general Eigen::MatrixXd into a linear operator
 */
template<typename Derived>
class EigenLinearOperator : public LinearOperator {
public:

  EigenLinearOperator(Eigen::MatrixBase<Derived> const& Ain) : LinearOperator(A.rows(), A.cols()), A(Ain){}

  virtual ~EigenLinearOperator(){};

  /** Apply the linear operator to a vector */
  virtual Eigen::MatrixXd Apply(Eigen::Ref<Eigen::MatrixXd> const& x) override {return A*x;};

  /** Apply the transpose of the linear operator to a vector. */
  virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<Eigen::MatrixXd> const& x) override {return A.transpose()*x};

  /** Fills in the reference \f$y\f$ with \f$y=Ax\f$ */
  virtual void Apply(Eigen::Ref<Eigen::MatrixXd> const& x, Eigen::Ref<Eigen::MatrixXd> y) override { y = A*x;};

  /** Fill in the reference \f$y\f$ with \f$y = A^Txf$ */
  virtual void ApplyTranspose(Eigen::Ref<Eigen::MatrixXd> const& x, Eigen::Ref<Eigen::MatrixXd> y) override {y = A.transpose()*x;};

protected:
  Eigen::MatrixBase<Derived> A;

};

template<typename Derived>
std::shared_ptr<LinearOperator> LinearOperator::Create(Eigen::MatrixBase<Derived> const& A)
{
    return std::make_shared<EigenLinearOperator>(A);
}


} // namespace Modeling
} // namespace MUQ



#endif
