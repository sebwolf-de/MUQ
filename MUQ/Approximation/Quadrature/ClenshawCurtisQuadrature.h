#ifndef CLENSHAWCURTISQUADRATURE_H
#define CLENSHAWCURTISQUADRATURE_H

#include "MUQ/Approximation/Quadrature/Quadrature.h"

namespace muq {
namespace Approximation {

  class ClenshawCurtisQuadrature : public Quadrature {
  public:

    ClenshawCurtisQuadrature();

    virtual ~ClenshawCurtisQuadrature() = default;

    virtual void Compute(unsigned int order) override;

  private:
    const double pi = 3.14159265358979323846;
  };

} // namespace muq
} // namespace Approximation


#endif
