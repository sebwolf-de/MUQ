#ifndef ODEDATA_H_
#define ODEDATA_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Modeling {
    class ODEData {
    public:

      ODEData(std::shared_ptr<WorkPiece> rhs, ref_vector<boost::any> const& inputs);

      ODEData(std::shared_ptr<WorkPiece> rhs, std::shared_ptr<WorkPiece> root, ref_vector<boost::any> const& inputs);

      std::shared_ptr<WorkPiece> rhs;

      std::shared_ptr<WorkPiece> root;

      ref_vector<boost::any> inputs;
      
    private:
    };
  } // namespace Modeling
} // namespace muq

#endif
