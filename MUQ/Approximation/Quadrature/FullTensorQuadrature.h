#ifndef FULLTENSORQUADRATURE_H_
#define FULLTENSORQUADRATURE_H_

#include <Eigen/Core>
#include <vector>

#include "MUQ/Approximation/Quadrature/Quadrature.h"


namespace muq {
  namespace Approximation {

/**
 * The simplest possible quadrature rule, a full tensor expansion.
 * Creates a the complete tensor product of the points of the given order
 * 1D quadrature rules for each dimension. The 1D quadrature order can
 * be isotropic or vary per dimension. Note that an isotropic order
 * does not mean isotropic number of points in each dimension, if different
 * 1D quadrature rules are used.
 *
 * Probably not the most useful class for actual analysis, because you have to know exactly
 * what order you want, but is an important building block for sparse
 * quadrature routines.
 */
class FullTensorQuadrature : public Quadrature {
public:

  FullTensorQuadrature(unsigned int                       dim,
                       std::shared_ptr<Quadrature> const& rules,
                       unsigned int                       order);

  FullTensorQuadrature(std::vector<std::shared_ptr<Quadrature>> const& rules,
                       Eigen::RowVectorXi                       const& orders);


  virtual ~FullTensorQuadrature() = default;

  virtual void Compute(unsigned int order) override;
  virtual void Compute(Eigen::RowVectorXi const& orders) override;

private:

  std::vector<std::shared_ptr<Quadrature>> rules;
};
}
}

#endif /* FULLTENSORQUADRATURE_H_ */
