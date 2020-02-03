#ifndef GAUSSNEWTONOPERATOR_H
#define GAUSSNEWTONOPERATOR_H

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/Distributions/GaussianBase.h"

#include <memory>

namespace muq
{
namespace Modeling
{


/** @class GaussNewtonOperator
 *  @ingroup LinearAlgebra
 *  @brief Creates a linear operator for the action of a Gauss-Newton Hessian approximation
           on a vector.
    @seealso HessianOperator
 */
class GaussNewtonOperator : public LinearOperator {
public:

  GaussNewtonOperator(std::shared_ptr<ModPiece>     const& forwardModelIn,
                      std::shared_ptr<GaussianBase> const& noiseModelIn,
                      std::vector<Eigen::VectorXd>  const& inputsIn,
                      unsigned int                         inWrt);

  virtual ~GaussNewtonOperator() = default;

  /** Apply the linear operator to a vector */
  virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

  /** Apply the transpose of the linear operator to a vector. */
  virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

protected:
  std::shared_ptr<ModPiece> forwardModel;
  std::shared_ptr<GaussianBase> noiseModel;

  const std::vector<Eigen::VectorXd> inputs;

  const unsigned int inWrt;

};

} // namespace Modeling
} // namespace MUQ



#endif
