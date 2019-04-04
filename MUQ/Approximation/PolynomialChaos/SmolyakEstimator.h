#ifndef SMOLYAKESTIMATOR_H
#define SMOLYAKESTIMATOR_H

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/Flann/FlannCache.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexSet.h"

#include <vector>
#include <boost/property_tree/ptree.hpp>

namespace muq {
namespace Approximation {

  template<typename EstimateType>
  class SmolyakEstimator {
  public:

    SmolyakEstimator(std::shared_ptr<muq::Modeling::ModPiece> const& modelIn);

    virtual ~SmolyakEstimator() = default;

    /** This is the main function to constructing static or adaptive Smolyak
        estimates.
    */
    virtual EstimateType Compute(std::shared_ptr<muq::Utilities::MultiIndexSet> const& fixedSet,
                                 boost::property_tree::ptree                           options = boost::property_tree::ptree());



  protected:

    /** Computes the locations where the model will need to be evaluated in order
        to construct a single tensor-product estimate.  For example, in the
        Smolyak quadrature setting, this function will return the quadrature points
        coming from the tensor product quadrature rule defined by the multiindex.

        This function works in tandem with the ComputeOneTerm function, which takes
        model evaluations and actually returns the tensor product estimate.  These
        two functions are split to allow the SmolyakEstimator to handle any job
        scheduling or caching that may be needed for parallel model evaluations.
    */
    virtual std::vector<Eigen::VectorXd> OneTermPoints(std::shared_ptr<muq::Utilities::MultiIndex> const& multi) = 0;

    /**
      This function works in tandem with the OneTermPoints function.  After the model
      has been evaluated at the points returned by OneTermPoints, this function will
      compute an estimate with the new model evaluations.  In the quadrature setting,
      the estimate will be computed with a weighted sum of the model evaluations stored
      in the modEvals vector.
    */
    virtual EstimateType ComputeOneTerm(std::shared_ptr<muq::Utilities::MultiIndex>                const& multi,
                                        std::vector<std::reference_wrapper<const Eigen::VectorXd>> const& modEvals) = 0;

    /** Should compute sum(smolyVals[i] * smolyWeights[i]) and return the result*/
    virtual EstimateType ComputeWeightedSum() const = 0;

    /** Computes the change in the smolyak weights caused by the addition of
        one term to the Smolyak MultiIndexSet.  This is to construct the c_k
        coefficients in equation 3.6 of Conrad and Marzouk's pseudospectral paper.
        @param[in] index The linear index of the new term in the "terms" set.
        @return A vector of length terms->Size() describing how the smolyWeights vector should be updated.
    */
    //void UpdateSmolyCoeffs(unsigned int const index);

    // /** Evaluates the model (possibly in parallel) at the specified points.  Also
    //     caches the results in the evalCache data structure to avoid any repeated
    //     model evaluations.
    // */
    // std::vector<Eigen::VectorXd> EvaluateModel(std::vector<Eigen::VectorXd> const& inputPts);

    /// The model used to construct the approximations
    std::shared_ptr<muq::Modeling::ModPiece> model;

    /// Multiindices defining each tensor product term in the Smolyak approximation
    std::shared_ptr<muq::Utilities::MultiIndexSet> terms;

    /// Contains each of the individual tensor product approximations
    std::vector<EstimateType> smolyVals;

    /// Holds the weights for the tensor product formulation.
    std::vector<double> smolyWeights;

    /// A cache of model evaluations
    muq::Modeling::DynamicKDTreeAdaptor<> pointCache;
    std::vector<Eigen::VectorXd> evalCache;


    int InCache(Eigen::VectorXd const& input) const;
    int AddToCache(Eigen::VectorXd const& newPt);
    int CacheSize() const{return pointCache.m_data.size();};

    const double cacheTol = 4e-15; // <- points are considered equal if they are closer than this

  }; // class SmolyakEstimator

} // namespace muq
} // namespace Approximation


#endif
