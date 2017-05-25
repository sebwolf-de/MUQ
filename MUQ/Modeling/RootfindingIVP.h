#ifndef ROOTFINDINGIVP_H_
#define ROOTFINDINGIVP_H_

#include <cvodes/cvodes.h> // prototypes for CVODE fcts. and consts.

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Modeling {
    class RootfindingIVP : public WorkPiece {
    public:

      RootfindingIVP();

      virtual ~RootfindingIVP();
      
    private:

      virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override;
      
    };
  } // namespace Modeling
} // namespace muq

#endif
