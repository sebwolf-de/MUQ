#ifndef FLANNCACHE_H_
#define FLANNCACHE_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Modeling {
    /// Create a cache of model evaluations (input/output pairs)
    /**
       Caches the input/output pairs for a muq::Modeling::WorkPiece that has one input and one output.
     */
    class FlannCache : public WorkPiece {
    public:

      /**
	 @param[in] function The function whose input/output pairs we want to cache
       */
      FlannCache(std::shared_ptr<WorkPiece> function);

      ~FlannCache();
      
    private:

      virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override;
      
    };
  } // namespace Utilities
} // namespace muq

#endif
