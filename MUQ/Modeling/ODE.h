#ifndef ODE_H_
#define ODE_H_

#include "MUQ/Modeling/ODEBase.h"

namespace muq {
  namespace Modeling {
    class ODE : public ODEBase {
    public:

      ODE(std::shared_ptr<WorkPiece> rhs);

      virtual ~ODE();

    private:

      virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override;
      
    };
  } // namespace Modeling
} // namespace muq

#endif
