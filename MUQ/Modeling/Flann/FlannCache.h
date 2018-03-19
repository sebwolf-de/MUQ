#ifndef FLANNCACHE_H_
#define FLANNCACHE_H_

//#include <flann/flann.hpp>
#include <nanoflann.hpp>

#include "MUQ/Utilities/LinearAlgebra/AnyAlgebra.h"

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Modeling {

    class FlannCache;


    template <class Distance = nanoflann::metric_L2, typename IndexType = size_t>
    struct DynamicKDTreeAdaptor
    {
        friend class FlannCache;

        typedef DynamicKDTreeAdaptor<Distance,IndexType> self_t;
        typedef typename Distance::template traits<double,self_t>::distance_t metric_t;
        typedef nanoflann::KDTreeSingleIndexDynamicAdaptor< metric_t,self_t,-1,IndexType>  index_t;

        std::shared_ptr<index_t> index; //! The kd-tree index for the user to call its methods as usual with any other FLANN index.
        std::vector<Eigen::VectorXd> m_data;

        /// Constructor: takes a const ref to the vector of vectors object with the data points
        DynamicKDTreeAdaptor(const int dim,
                             const int leaf_max_size = 10)
        {
          index = std::make_shared<index_t>(dim, *this /* adaptor */, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size ) );
        }

        DynamicKDTreeAdaptor(std::vector<Eigen::VectorXd> const& pts,
                             const int leaf_max_size = 10) : m_data(pts)
        {
          index = std::make_shared<index_t>(pts.at(0).size(), *this /* adaptor */, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size ) );
          index->addPoints(0, pts.size());
        }

        inline void add(Eigen::VectorXd const& newPt)
        {
          m_data.push_back(newPt);
          index->addPoints(m_data.size()-1, m_data.size()-1);
        }

        /** Query for the \a num_closest closest points to a given point (entered as query_point[0:dim-1]).
          *  Note that this is a short-cut method for index->findNeighbors().
          *  The user can also call index->... methods as desired.
          * \note nChecks_IGNORED is ignored but kept for compatibility with the original FLANN interface.
          */
        inline std::pair<std::vector<IndexType>, std::vector<double>> query(Eigen::VectorXd const& query_point,
                                                                            const size_t num_closest,
                                                                            const int nChecks_IGNORED = 10) const
        {
          std::vector<IndexType> out_indices(num_closest);
          std::vector<double> out_distances_sq(num_closest);

          nanoflann::KNNResultSet<double,IndexType> resultSet(num_closest);
          resultSet.init(&out_indices[0], &out_distances_sq[0]);
          index->findNeighbors(resultSet, query_point.data(), nanoflann::SearchParams());

          return std::make_pair(out_indices, out_distances_sq);
        }

        /** @name Interface expected by KDTreeSingleIndexAdaptor
          * @{ */

        const self_t & derived() const {
          return *this;
        }
        self_t & derived()       {
          return *this;
        }

        // Must return the number of data points
        inline size_t kdtree_get_point_count() const {
          return m_data.size();
        }

        // Returns the dim'th component of the idx'th point in the class:
        inline double kdtree_get_pt(const size_t idx, int dim) const {
          return m_data[idx][dim];
        }

        // Optional bounding-box computation: return false to default to a standard bbox computation loop.
        //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
        //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
        template <class BBOX>
        bool kdtree_get_bbox(BBOX & /*bb*/) const {
          return false;
        }

    }; // end of DynamicKDTreeAdaptor


    /// Create a cache of model evaluations (input/output pairs)
    /**
       Caches the input/output pairs for a muq::Modeling::WorkPiece that has one input and one output.
     */
    class FlannCache : public WorkPiece {
    public:

      /**
	     @param[in] function The function whose input/output pairs we want to cache
       */
      FlannCache(std::shared_ptr<WorkPiece> function,
                 int                        inputDim);

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
      void Add(Eigen::VectorXd const& input);

      /// Remove point from the cache
      /**
	 @param[in] input The entry we would like to remove from the cache
       */
      void Remove(boost::any const& input);

      /// Find the \f$k\f$ nearest neighbors
      /**
	 @param[in] point The point whose nearest neighbors we want to find
	 @param[in] k We want to find this many nearest neighbors
	 @param[out] neighbors A vector of the \fk\f$ nearest neighbors
	 @param[out] result The output corresponding to the \f$k\f$ nearest neighbors
       */
      void NearestNeighbors(boost::any const& point,
                            unsigned int const k,
                            std::vector<Eigen::VectorXd>& neighbors,
                            std::vector<Eigen::VectorXd>& result) const;

      /// Get the size of the cache
      /**
	     \return The size of the cache
       */
      unsigned int Size() const;

    private:

      // The vector of previous results
      std::vector<Eigen::VectorXd> outputCache;

      virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override;

      /// The function whose input/outputs we are caching
      std::shared_ptr<WorkPiece> function;

      /// The nearest neighbor index, used to perform searches
      std::shared_ptr<DynamicKDTreeAdaptor<>> kdTree;

    };
  } // namespace Utilities
} // namespace muq

#endif
