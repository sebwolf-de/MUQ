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

  unsigned numQuad = 200;
  quadPts = Eigen::VectorXd::LinSpaced(numQuad+1,0,1).head(numQuad);
  quadWeights = (1.0/numQuad)*Eigen::VectorXd::Ones(numQuad);
}


void MonotoneExpansion::EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs)
{

  bool updateCoeffs = inputs.size()>1;

  Eigen::VectorXd coeffs;
  if(updateCoeffs){
    coeffs = boost::any_cast<Eigen::VectorXd>(inputs.at(1).get());
  }

  Eigen::VectorXd output = Eigen::VectorXd::Zero(monotoneParts.size());

  Eigen::VectorXd const& x = *boost::any_cast<Eigen::VectorXd>(&inputs.at(0).get());

  Eigen::MatrixXd tempCoeffs;
  unsigned currCoeff = 0;

  // Fill in the general parts
  for(int i=0; i<generalParts.size(); ++i){

    if(updateCoeffs){
      tempCoeffs = coeffs.segment(currCoeff, generalParts.at(i)->NumTerms()).transpose();
      output(i+1) += boost::any_cast<Eigen::VectorXd>(generalParts.at(i)->Evaluate(x.head(i+1).eval(), tempCoeffs).at(0))(0);
    }else{
      output(i+1) += boost::any_cast<Eigen::VectorXd>(generalParts.at(i)->Evaluate(x.head(i+1).eval()).at(0))(0);
    }

    currCoeff += generalParts.at(i)->NumTerms();
  }

  // Now add the monotone part
  Eigen::VectorXd quadEvals(quadPts.size());
  for(int i=0; i<monotoneParts.size(); ++i){


    // Loop over quadrature points
    Eigen::VectorXd evalPt = x.head(i+1);
    evalPt(i) = x(i) * quadPts(0);

    double polyEval;
    if(updateCoeffs){
      tempCoeffs = coeffs.segment(currCoeff, monotoneParts.at(i)->NumTerms()).transpose();
      polyEval = boost::any_cast<Eigen::VectorXd>(monotoneParts.at(i)->Evaluate(evalPt, tempCoeffs).at(0))(0);
    }else{
      polyEval = boost::any_cast<Eigen::VectorXd>(monotoneParts.at(i)->Evaluate(evalPt).at(0))(0);
    }
    currCoeff += monotoneParts.at(i)->NumTerms();

    quadEvals(0) = polyEval*polyEval;
    for(int k=1; k<quadPts.size(); ++k){
      evalPt(i) = x(i) * quadPts(k);
      polyEval = boost::any_cast<Eigen::VectorXd>(monotoneParts.at(i)->Evaluate(evalPt).at(0))(0);
      quadEvals(k) = polyEval*polyEval;
    }
    output(i) += quadEvals.dot(quadWeights)*x(i);
  }

  outputs.resize(1);
  outputs.at(0) = output;
}

void MonotoneExpansion::JacobianImpl(unsigned int const                           wrtIn,
                                     unsigned int const                           wrtOut,
                                     muq::Modeling::ref_vector<boost::any> const& inputs)
{

}
