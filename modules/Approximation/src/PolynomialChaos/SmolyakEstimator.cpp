#include "MUQ/Approximation/PolynomialChaos/SmolyakEstimator.h"

#include "MUQ/Approximation/Quadrature/SmolyakQuadrature.h"
#include "MUQ/Approximation/PolynomialChaos/PolynomialChaosExpansion.h"

using namespace muq::Modeling;
using namespace muq::Utilities;
using namespace muq::Approximation;

template<typename EstimateType>
SmolyakEstimator<EstimateType>::SmolyakEstimator(std::shared_ptr<muq::Modeling::ModPiece> const& modelIn)
  : model(modelIn),
    terms(std::make_shared<MultiIndexSet>(modelIn->inputSizes(0))),
    pointCache(modelIn->inputSizes(0))
{
  // make sure the model has a single input and output
  assert(model->inputSizes.size()==1);
  assert(model->outputSizes.size()==1);
}


template<typename EstimateType>
EstimateType SmolyakEstimator<EstimateType>::Compute(std::shared_ptr<muq::Utilities::MultiIndexSet> const& fixedSet,
                                                    boost::property_tree::ptree                           options)
{
  /* For each term in the fixed set we want to construct a tensor product approximation.
     The i^th tensor product approximation will require model evaluations at N_i
     points.   The allNewPts vector will store all of the points that need to be evaluated
     and stored in our cache.  The evalInds vector will keep track of which points
     were needed for each tensor product approximation by storing the index of
     the cached points.
  */
  std::vector<std::vector<unsigned int>> evalInds(fixedSet->Size());
  std::vector<Eigen::VectorXd> allNewPts; // new pts to evaluate

  for(unsigned int i=0; i<fixedSet->Size(); ++i) {

    std::vector<Eigen::VectorXd> pts = OneTermPoints(fixedSet->IndexToMulti(i));

    // Figure out if we've already evaluated a point, or if we need to
    evalInds.at(i).resize(pts.size());
    for(unsigned int k=0; k<pts.size(); ++k){

      // Check if the point is already in the cache
      int cacheId = InCache(pts.at(k));
      if(cacheId<0){
        int cacheInd = AddToCache(pts.at(k));
        evalInds.at(i).at(k) = cacheInd;
        allNewPts.push_back(pts.at(k));
      }else{
        evalInds.at(i).at(k) = cacheId;
      }
    }
  }

  /* The allNewPts vector now contains all the new points we need to evaluate and add to our cache
     In the future, this for loop could easily be parallelized
  */
  for(unsigned int ptInd=0; ptInd<allNewPts.size(); ++ptInd)
    evalCache.push_back( model->Evaluate(allNewPts.at(ptInd)).at(0) );

  // Resize the Smolyak weights to allow for the new terms
  smolyWeights.resize(smolyWeights.size() + fixedSet->Size(), 0.0);

  // Now, pull out references to each output vector and compute the tensor product estimates
  for(unsigned int i=0; i<fixedSet->Size(); ++i) {

    // Copy references of the model output to a vector
    std::vector<std::reference_wrapper<const Eigen::VectorXd>> evals;
    for(unsigned int ptInd=0; ptInd<evalInds.at(i).size(); ++ptInd){
      evals.push_back( evalCache.at(evalInds.at(i).at(ptInd)) );
    }

    // Compute the tensor product estimate using the model output
    int newInd = terms->AddActive(fixedSet->at(i));

    smolyVals.push_back( ComputeOneTerm(fixedSet->IndexToMulti(i), evals) );
    assert(newInd==smolyVals.size()-1);

    // update the smolyak coefficients now that we've added this term
    //UpdateSmolyCoeffs(newInd);

  }

  auto tempWeights = SmolyakQuadrature::ComputeWeights(terms);
  for(unsigned int i=0; i<tempWeights.size(); ++i)
    smolyWeights.at(i) = tempWeights(i);

  // Make sure nothing funky happened
  assert(smolyVals.size() == smolyWeights.size());

  // We've done all the work, just return a weighted sum of the tensor product approximations
  return ComputeWeightedSum();
}


// template<typename EstimateType>
// void SmolyakEstimator<EstimateType>::UpdateSmolyCoeffs(unsigned int const index)
// {
//   // Get backward neighbors
//   std::vector<unsigned int> neighInds = terms->GetBackwardNeighbors(index);
//
//   Eigen::RowVectorXi baseMulti = terms->IndexToMulti(index)->GetVector();
//   smolyWeights.at(index) += 1;
//
//   for(unsigned int neighInd : neighInds) {
//     int parity = (baseMulti - terms->IndexToMulti(neighInd)->GetVector()).sum();
//     smolyWeights.at(neighInd) += ( (parity % 2 == 0) ? 1 : -1);
//   }
// }

template<typename EstimateType>
int SmolyakEstimator<EstimateType>::AddToCache(Eigen::VectorXd const& newPt)
{
  int cacheId = InCache(newPt);
  if(cacheId<0){
    pointCache.add(newPt);
    return CacheSize()-1;
  }else{
    return cacheId;
  }
}

template<typename EstimateType>
int SmolyakEstimator<EstimateType>::SmolyakEstimator::InCache(Eigen::VectorXd const& input) const
{
  if(CacheSize()>0){
    std::vector<size_t> indices;
    std::vector<double> squaredDists;
    std::tie(indices, squaredDists) = pointCache.query(input, 1);

    if(squaredDists.at(0)<cacheTol){
      return indices.at(0);
    }else{
      return -1;
    }

  }else{
    return -1;
  }
}

namespace muq{
namespace Approximation{
  template class SmolyakEstimator<Eigen::VectorXd>;
  template class SmolyakEstimator<std::shared_ptr<PolynomialChaosExpansion>>;
}
}
