#ifndef LOCALREGRESSION_H_
#define LOCALREGRESSION_H_

#include "MUQ/Modeling/Flann/FlannCache.h"

namespace muq {
  namespace Approximation {
    class LocalRegression : public muq::Modeling::WorkPiece {
    public:

      /**
	 @param[in] function The function we wish to approximate with a local polynomial
       */
      LocalRegression(std::shared_ptr<muq::Modeling::WorkPiece> function);

      ~LocalRegression();

      /// Add some points 

    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      /// A cache containing previous model evaluations
      std::shared_ptr<muq::Modeling::FlannCache> cache;
      
    };
  } // namespace Approximation
} // namespace muq

#endif
