#include "MUQ/Approximation/Quadrature/AdaptiveSmolyakQuadrature.h"

using namespace muq::Modeling;
using namespace muq::Approximation;
using namespace muq::Utilities;


AdaptiveSmolyakQuadrature::AdaptiveSmolyakQuadrature(std::shared_ptr<muq::Modeling::ModPiece> const& modelIn,
                                                     std::vector<std::shared_ptr<Quadrature>> const& quad1d)
            : SmolyakEstimator<Eigen::VectorXd>(modelIn),
              tensQuad(quad1d)

{
  assert(modelIn->inputSizes.size()==1);
  assert(modelIn->inputSizes(0)==quad1d.size());
}


std::vector<Eigen::VectorXd> AdaptiveSmolyakQuadrature::OneTermPoints(std::shared_ptr<MultiIndex> const& multi)
{
  tensQuad.Compute(multi->GetVector());
  cachedMulti = multi;

  std::vector<Eigen::VectorXd> results(tensQuad.Points().cols());
  for(unsigned int i=0; i<tensQuad.Points().cols(); ++i)
    results.at(i) = tensQuad.Points().col(i);

  return results;
}


Eigen::VectorXd AdaptiveSmolyakQuadrature::ComputeOneTerm(std::shared_ptr<MultiIndex>                                const& multi,
                                                          std::vector<std::reference_wrapper<const Eigen::VectorXd>> const& modEvals)
{
  if((*cachedMulti) != (*multi)){
    tensQuad.Compute(multi->GetVector());
    cachedMulti = multi;
  }

  Eigen::VectorXd const& wts = tensQuad.Weights();
  assert(wts.size()==modEvals.size());

  // Compute the quadrature
  Eigen::VectorXd output = modEvals.at(0).get() * wts(0);
  for(unsigned int i=1; i<wts.size(); ++i)
    output += modEvals.at(i).get() * wts(i);

  return output;
}


Eigen::VectorXd AdaptiveSmolyakQuadrature::ComputeWeightedSum() const
{
  unsigned int firstGoodInd;
  for(unsigned int i=0; i<smolyWeights.size(); ++i){
    if(std::abs(smolyWeights.at(i))>1e-10){
      firstGoodInd = i;
      break;
    }
  }
  Eigen::VectorXd output = smolyWeights.at(firstGoodInd)*smolyVals.at(firstGoodInd);
  for(unsigned int i=firstGoodInd+1; i<smolyVals.size(); ++i){
    if(std::abs(smolyWeights.at(i))>1e-10)
      output += smolyWeights.at(i) * smolyVals.at(i);
  }

  return output;
}
