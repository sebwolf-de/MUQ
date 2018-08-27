#ifndef LOCALREGRESSION_H_
#define LOCALREGRESSION_H_

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/Flann/FlannCache.h"

#include "MUQ/Approximation/Regression/Regression.h"

namespace muq {
  namespace Approximation {
    class LocalRegression : public muq::Modeling::ModPiece {
    public:

      /**
	 @param[in] function The function we wish to approximate with a local polynomial
       */
      LocalRegression(std::shared_ptr<muq::Modeling::ModPiece> function, boost::property_tree::ptree& pt);

      ~LocalRegression();

      /// Add some points to the cache
      /**
	 @param[in] inputs Points to add to the cache
       */
      void Add(std::vector<Eigen::VectorXd> const& inputs) const;

      /// Add a single point to the cache
      /**
	 @param[in] input Point to add to the cache
	 \return The function result at the new point
       */
      Eigen::VectorXd Add(Eigen::VectorXd const& input) const;

      /// Get the total size of the cache
      /**
	 \return The cache size
       */
      unsigned int CacheSize() const;

      /// Get the poisedness constant
      /**
	 Get the poisedness constant associated with the nearest neighbors of the input point.
	 @param[in] input The input point
	 \return first: the point where the Lagrange polynomials are maximized, second: the poisedness constant
       */
      std::pair<Eigen::VectorXd, double> PoisednessConstant(Eigen::VectorXd const& input) const;

      /// Get the poisedness constant
      /**
	 Get the poisedness constant associated with the nearest neighbors of the input point.
	 @param[in] input The input point
	 @param[in] neighbors The nearest neighbors
	 \return first: the point where the Lagrange polynomials are maximized, second: the poisedness constant
       */
      std::pair<Eigen::VectorXd, double> PoisednessConstant(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd> const& neighbors) const;

      /// Get the error indicator
      /**
	 Get the error indicator \f$\Lambda \sqrt{k} \Delta^{p+1}\f$
	 @param[in] input The input point
	 \return first: the point where the Lagrange polynomials are maximized, second: the error indicator
       */
      std::pair<Eigen::VectorXd, double> ErrorIndicator(Eigen::VectorXd const& input) const;

      /// Get the error indicator
      /**
	 Get the error indicator \f$\Lambda \sqrt{k} \Delta^{p+1}\f$
	 @param[in] input The input point
	 @param[in] neighbors The nearest neighbors
	 \return first: the point where the Lagrange polynomials are maximized, second: the error indicator
       */
      std::pair<Eigen::VectorXd, double> ErrorIndicator(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd> const& neighbors) const;

      /// Get the number of nearest neighbors
      /**
	 @param[in] input We want the \f$k\f$ nearest neighbors to this point
	 @param[out] neighbors The \f$k\f$ nearest neighbors
       */
      void NearestNeighbors(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd>& neighbors) const;

      /// Get the number of nearest neighbors (with result)
      /**
	 @param[in] input We want the \f$k\f$ nearest neighbors to this point
	 @param[out] neighbors The \f$k\f$ nearest neighbors
	 @param[out] results The corresponding output of the function
       */
      void NearestNeighbors(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd>& neighbors, std::vector<Eigen::VectorXd>& result) const;

      Eigen::VectorXd EvaluateRegressor(Eigen::VectorXd const& input, std::vector<Eigen::VectorXd> const& neighbors, std::vector<Eigen::VectorXd> const& result) const;

      /// The number of nearest neighbors to use by the regressor
      const unsigned int kn;

    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override;

      /// Fit the regression to the nearest neighbors
      void FitRegression(Eigen::VectorXd const& input) const;

      /// A cache containing previous model evaluations
      std::shared_ptr<muq::Modeling::FlannCache> cache;

      /// A regressor
      std::shared_ptr<Regression> reg;      
    };
  } // namespace Approximation
} // namespace muq

#endif
