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
  std::vector<Eigen::VectorXd> allNewPts; // new pts to evaluate

  for(unsigned int i=0; i<fixedSet->Size(); ++i) {
    evalInds.push_back(std::vector<unsigned int>());
    terms->AddActive(fixedSet->at(i));
    smolyVals.push_back( EstimateType());
    isComputed.push_back(false);

    assert(evalInds.size()==terms->Size());

    std::vector<Eigen::VectorXd> pts = OneTermPoints(fixedSet->IndexToMulti(i));

    // Figure out if we've already evaluated a point, or if we need to
    evalInds.at(i).resize(pts.size());
    for(unsigned int k=0; k<pts.size(); ++k){

      // Check if the point is already in the cache
      int cacheId = InCache(pts.at(k));
      if(cacheId<0){
        int cacheInd = AddToCache(pts.at(k));
        evalInds.at(i).at(k) = cacheInd;
        evalCache.push_back( Eigen::VectorXd() );

      }else{
        evalInds.at(i).at(k) = cacheId;
      }
    }
  }

  // Now, compute the Smolyak coefficients for the new multiindex set
  auto tempWeights = SmolyakQuadrature::ComputeWeights(terms);
  smolyWeights.resize(tempWeights.size());
  for(unsigned int i=0; i<tempWeights.size(); ++i)
    smolyWeights.at(i) = tempWeights(i);

  // Figure out what points we need to evaluate
  std::set<unsigned int> ptsToEval;
  for(unsigned int termInd=0; termInd<evalInds.size(); ++termInd) {

    // If the smolyak weight is zero, don't bother computing any needed points
    if( std::abs(smolyWeights.at(termInd)) > 10*std::numeric_limits<double>::epsilon()) {

      for(unsigned int i=0; i<evalInds.at(termInd).size(); ++i) {
        unsigned int ptInd = evalInds.at(termInd).at(i);
        // If the size of the output is zero, we haven't evaluated the model at this point yet
        if(evalCache.at(ptInd).size()==0)
          ptsToEval.insert(ptInd);
      }
    }
  }

  // Evaluate all the terms we need to
  for(auto& ptInd : ptsToEval)
    evalCache.at(ptInd) = model->Evaluate( GetFromCache(ptInd) ).at(0);

  // Now, compute any new tensor product estimates that are necessary
  for(unsigned int i=0; i<terms->Size(); ++i) {

    // If we haven't built this estimate before, do it now
    if((!isComputed.at(i))&&( std::abs(smolyWeights.at(i))>1e-10)) {

      // Copy references of the model output to a vector
      std::vector<std::reference_wrapper<const Eigen::VectorXd>> evals;
      for(unsigned int ptInd=0; ptInd<evalInds.at(i).size(); ++ptInd){
        evals.push_back( evalCache.at(evalInds.at(i).at(ptInd) ) );
      }
      smolyVals.at(i) = ComputeOneTerm(fixedSet->IndexToMulti(i), evals);
    }

  }

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
