#ifndef LOCALREGRESSION_H_
#define LOCALREGRESSION_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Approximation {
    class LocalRegression : public muq::Modeling::WorkPiece {
    public:

      /**
	 @param[in] function The function we wish to approximate with a local polynomial
       */
      LocalRegression(std::shared_ptr<muq::Modeling::WorkPiece> function);

      ~LocalRegression();
      
    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;
      
    };
  } // namespace Approximation
} // namespace muq

#endif
