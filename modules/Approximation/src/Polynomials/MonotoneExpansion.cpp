#include "MUQ/Approximation/Polynomials/MonotoneExpansion.h"


using namespace muq::Approximation;
using namespace muq::Utilities;

MonotoneExpansion::MonotoneExpansion(std::shared_ptr<BasisExpansion> monotonePartsIn) : MonotoneExpansion({},{monotonePartsIn})
{
}

MonotoneExpansion::MonotoneExpansion(std::vector<std::shared_ptr<BasisExpansion>> const& generalPartsIn,
                                     std::vector<std::shared_ptr<BasisExpansion>> const& monotonePartsIn) : generalParts(generalPartsIn),
                                                                                                            monotoneParts(monotonePartsIn)
{
  assert(generalPartsIn.size() == monotonePartsIn.size()-1);

  unsigned numQuad = 20;
  quadPts = Eigen::VectorXd::LinSpaced(0,1,numQuad+1).head(numQuad);
  quadWeights = (1.0/numQuad)*Eigen::VectorXd::Ones(numQuad);
}


void MonotoneExpansion::EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs)
{

  Eigen::VectorXd output = Eigen::VectorXd::Zero(monotoneParts.size());

  Eigen::VectorXd const& x = *boost::any_cast<Eigen::VectorXd>(&inputs.at(0).get());

  Eigen::VectorXd tempCoeffs;
  unsigned currCoeff = 0;

  // Fill in the general parts
  for(int i=0; i<generalParts.size(); ++i){
    tempCoeffs = coeffs.block(currCoeff, 0, 1, generalParts.at(i)->NumTerms());
    output(i+1) += generalParts.at(i)->Evaluate(x.head(i+1).eval(), tempCoeffs);

    currCoeff += generalParts.at(i)->NumTerms();
  }

  // Now add the monotone part
  Eigen::VectorXd quadEvals(quadPts.size());
  for(int i=0; i<monotoneParts.size(); ++i){

    // Loop over quadrature point
    Eigen::VectorXd evalPt = x.head(i+1);
    for(int k=0; k<quadPts.size(); ++k){
      evalPt(i) = x(i) * quadPts;
      tempCoeffs = coeffs.block(currCoeff, 0, 1, monotoneParts.at(i)->NumTerms());
      quadEvals(k) = monotoneParts->Evaluate(evalPt, tempCoeffs);

      currCoeff += monotoneParts.at(i)->NumTerms();
    }

    output(i) += x(i)*quadEvals.dot(quadWeights);
  }


  outputs.resize(1);
  outputs.at(0) = output;
}

void MonotoneExpansion::JacobianImpl(unsigned int const                           wrtIn,
                                     unsigned int const                           wrtOut,
                                     muq::Modeling::ref_vector<boost::any> const& inputs)
{

}
