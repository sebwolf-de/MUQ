#ifndef LOCALREGRESSION_H_
#define LOCALREGRESSION_H_

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Modeling/Flann/FlannCache.h"

#include "MUQ/Approximation/Regression/Regression.h"

namespace muq {
  namespace Approximation {
    class LocalRegression : public muq::Modeling::WorkPiece {
    public:

      /**
	 @param[in] function The function we wish to approximate with a local polynomial
       */
      LocalRegression(std::shared_ptr<muq::Modeling::WorkPiece> function, boost::property_tree::ptree const& pt);

      ~LocalRegression();

      /// Add some points to the cache
      template<typename intype>
	inline void Add(std::vector<intype> const& inputs) const {
	assert(cache);
	cache->Add<intype>(inputs);
      }

      /// Get the total size of the cache
      /**
	 \return The cache size
       */
      unsigned int CacheSize() const;

    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      /// Fit the regression to the nearest neighbors
      void FitRegression(boost::any const& input) const;

      /// A cache containing previous model evaluations
      std::shared_ptr<muq::Modeling::FlannCache> cache;

      /// A regressor
      std::shared_ptr<Regression> reg;

      /// The number of nearest neighbors to use by the regressor
      const unsigned int kn;
      
    };
  } // namespace Approximation
} // namespace muq

#endif
