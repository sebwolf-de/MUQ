#ifndef FLANNCACHE_H_
#define FLANNCACHE_H_

#include <flann/flann.hpp>

#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Modeling/AnyAlgebra.h"

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

      /// Determine if an entry is in the cache
      /**
	 @param[in] input Check if this input is in the cache
	 \return -1: the input vector does not have an entry in the cache, otherise: return the index of the input in the cache
       */
      int InCache(boost::any const& input) const;

      /// Add new points to the cache
      template<typename intype>
	inline void Add(std::vector<intype> const& inputs) {
	for( auto it : inputs ) {
	  // add the point if is not already there
	  if( InCache(it)<0 ) { Add(it); }

	  // make sure it got added
	  assert(InCache(it)>=0);
	}
      }

      /// Add a new point to the cache
      /**
	 @param[in] input The entry we would like to add to the cache (if it is not there already)
       */
      void Add(boost::any const& input);

      /// Remove point from the cache
      /**
	 @param[in] input The entry we would like to remove from the cache
       */
      void Remove(boost::any const& input);

      /// Find the \f$k\f$ nearest neighbors
      /**
	 @param[in] point The point whose nearest neighbors we want to find
	 @param[in] k We want to find this many nearest neighbors	 
       */
      void NearestNeighbors(boost::any const& point, unsigned int const k) const;

      /// Get the size of the cache
      /**
	 \return The size of the cache
       */
      unsigned int Size() const;
      
    private:

      /// An entry in the cache
      /**
	 Contains: input, output, derivative information
       */
      struct Entry {
	Entry();

	/// The output for this entry
	flann::Matrix<double> output;
      };

      /// Copy any vector into a flann type vector
      /**
	 @param[in] vec The vector we wish to convert to a flann type
	 @param[out] fvec The copied result
       */
      void DeepVectorCopy(boost::any const& vec, flann::Matrix<double>& fvec) const;
      
      virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override;

      /// The function whose input/outputs we are caching
      std::shared_ptr<WorkPiece> function;

      /// An algebra for matrix/vector manipulation
      std::shared_ptr<AnyAlgebra> algebra;

      /// The cache---maps unique ids to the FlannCache::Entry
      std::map<unsigned int, std::shared_ptr<Entry> > cache;

      /// The nearest neighbor index, used to perform searches
      std::shared_ptr<flann::Index<flann::L2<double> > > nnIndex;

      // the id of the next point that will be added
      unsigned int nextID = 0; 
      
    };
  } // namespace Utilities
} // namespace muq

#endif
