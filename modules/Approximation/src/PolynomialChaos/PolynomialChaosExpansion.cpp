#include "MUQ/Approximation/PolynomialChaos/PolynomialChaosExpansion.h"

#include <fstream>
#include <iostream>

#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

using namespace muq::Utilities;
using namespace muq::Approximation;

PolynomialChaosExpansion::PolynomialChaosExpansion(std::shared_ptr<OrthogonalPolynomial>          const& basisCompsIn,
                                                   std::shared_ptr<muq::Utilities::MultiIndexSet>        multisIn,
                                                   Eigen::MatrixXd                                const& coeffsIn) : PolynomialChaosExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>>(multisIn->GetMultiLength(), basisCompsIn),
                                                                                                                                              multisIn,
                                                                                                                                              coeffsIn)
{}

PolynomialChaosExpansion::PolynomialChaosExpansion(std::shared_ptr<OrthogonalPolynomial>          const& basisCompsIn,
                                                   std::shared_ptr<muq::Utilities::MultiIndexSet>        multisIn,
                                                   unsigned int                                          outputDim) : PolynomialChaosExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>>(multisIn->GetMultiLength(), basisCompsIn),
                                                                                                                                               multisIn,
                                                                                                                                               outputDim)
{}

PolynomialChaosExpansion::PolynomialChaosExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                                                   std::shared_ptr<muq::Utilities::MultiIndexSet>          multisIn,
                                                   Eigen::MatrixXd                                  const& coeffsIn) : BasisExpansion(basisCompsIn, multisIn, coeffsIn)
{
}

PolynomialChaosExpansion::PolynomialChaosExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                                                   std::shared_ptr<muq::Utilities::MultiIndexSet>            multisIn,
                                                   unsigned int                                              outputDim) : BasisExpansion(basisCompsIn, multisIn, Eigen::MatrixXd::Zero(outputDim, multisIn->GetMultiLength()))
{}



Eigen::VectorXd PolynomialChaosExpansion::GetNormalizationVec() const{

  Eigen::VectorXd result = Eigen::VectorXd::Zero(multis->Size());

  //compute each one
  for (unsigned int i = 0; i < multis->Size(); i++) {
    double norm = 1.0; //start the normalization at 1

    //loop over the dimensions and multiply the normalizations for each component polynomial
    std::shared_ptr<MultiIndex> multi = multis->IndexToMulti(i);
    for(auto it=multi->GetNzBegin(); it!=multi->GetNzEnd(); ++it)
       norm *= std::dynamic_pointer_cast<OrthogonalPolynomial>(basisComps.at(it->first))->Normalization(it->second);

    result(i) = std::sqrt(norm);
  }
  return result;
};


Eigen::VectorXd PolynomialChaosExpansion::Variance() const
{
  Eigen::VectorXd normalVec = GetNormalizationVec();

  //if there's only the one constant term, the PCE has variance zero
  if (normalVec.rows() <= 1) {
    return Eigen::VectorXd::Zero(coeffs.rows());
  }

  Eigen::VectorXd squareNorm = normalVec.tail(normalVec.rows() - 1).array().square();

  //for each output, the variance is the dot product of the squared coeffs with the squared norms of the PCE terms.
  //Thus, grab all but the leftmost column.
  //Have to normalize by what the constant integrates to.

  Eigen::VectorXd dimMeasures(inputSizes(0));
  for (unsigned int i = 0; i < inputSizes(0); ++i)
    dimMeasures(i) = std::dynamic_pointer_cast<OrthogonalPolynomial>(basisComps.at(i))->Normalization(0);

  //Since the variance is an expectation, we must normalize if the polynomials aren't quite set up
  //to integrate to one.
  return coeffs.rightCols(coeffs.cols() - 1).array().square().matrix() * squareNorm / dimMeasures.prod();
}

Eigen::MatrixXd PolynomialChaosExpansion::Covariance() const
{
  Eigen::VectorXd normalVec = GetNormalizationVec();

  //if there's only the one constant term, the PCE has variance zero
  if (normalVec.rows() <= 1) {
    return Eigen::MatrixXd::Zero(coeffs.rows(),coeffs.rows());
  }

  Eigen::VectorXd squareNorm = normalVec.tail(normalVec.rows() - 1).array().square();

  //for each output, the variance is the dot product of the squared coeffs with the squared norms of the PCE terms.
  //Thus, grab all but the leftmost column.
  //Have to normalize by what the constant integrates to.

  Eigen::VectorXd dimMeasures(inputSizes(0));
  for (unsigned int i = 0; i < inputSizes(0); ++i)
    dimMeasures(i) = std::dynamic_pointer_cast<OrthogonalPolynomial>(basisComps.at(i))->Normalization(0);

  Eigen::VectorXd quadVec = squareNorm / dimMeasures.prod();

  // fill in upper portion of covariance matrix
  Eigen::MatrixXd cov(outputSizes(0), outputSizes(0));
  int npolys = coeffs.cols();
  for (int j = 0; j < outputSizes(0); ++j) {   // column index
    for (int i = j; i <= outputSizes(0); ++i) { // row index
      cov(i,j) = coeffs.row(j).tail(npolys - 1) * quadVec.asDiagonal() * coeffs.row(i).tail(npolys - 1).transpose();
      cov(j,i) = cov(i,j);
    }
  }

  // lower to upper triangular copy
  for (int j = 1; j < outputSizes(0); ++j) { // column index
    for (int i = 0; i < j; ++i) {        // row index
      cov(i, j) = cov(j, i);
    }
  }

  return cov;
}

Eigen::VectorXd PolynomialChaosExpansion::Mean() const
{
  return coeffs.col(0);
}

// Eigen::MatrixXd PolynomialChaosExpansion::ComputeVarianceJacobian() const
// {
//   VectorXd normalVec = GetNormalizationVec();
//
//   //if there's only the one constant term, the PCE has variance zero
//   if (normalVec.rows() <= 1) {
//     return MatrixXd::Zero(coeffs.rows(), coeffs.cols());
//   }
//
//   unsigned int npolys = coeffs.cols();
//
//   VectorXd squareNorm = normalVec.tail(normalVec.rows() - 1).array().square();
//
//   VectorXd dimMeasures(inputSizes(0));
//   for (unsigned int i = 0; i < inputSizes(0); ++i)
//     dimMeasures(i) = polys[i]->GetNormalization(0);
//
//   Eigen::MatrixXd Jacobian = Eigen::MatrixXd::Zero(coeffs.rows(), npolys);
//
//   Eigen::VectorXd scales = squareNorm / dimMeasures.prod();
//   for (unsigned int i = 1; i < npolys; ++i) {
//     Jacobian.col(i) = 2.0 * coeffs.col(i) * scales(i - 1);
//   }
//
//   //Since the variance is an expectation, we must normalize if the polynomials aren't quite set up
//   //to integrate to one.
//   return Jacobian;
// }
//
// Eigen::VectorXd PolynomialChaosExpansion::ComputeVarianceGradient(int outInd) const
// {
//   VectorXd normalVec = GetNormalizationVec();
//
//   //if there's only the one constant term, the PCE has variance zero
//   if (normalVec.rows() <= 1) {
//     return VectorXd::Zero(1);
//   }
//
//   VectorXd squareNorm = normalVec.tail(normalVec.rows() - 1).array().square();
//
//   VectorXd dimMeasures(inputSizes(0));
//   for (unsigned int i = 0; i < inputSizes(0); ++i) {
//     dimMeasures(i) = polys[i]->GetNormalization(0);
//   }
//
//   Eigen::VectorXd Gradient = Eigen::VectorXd::Zero(coeffs.cols());
//
//   Eigen::VectorXd scales = squareNorm / dimMeasures.prod();
//   for (unsigned int i = 1; i < coeffs.cols(); ++i) {
//     Gradient(i) = coeffs(outInd, i) * scales(i - 1);
//   }
//
//   //Since the variance is an expectation, we must normalize if the polynomials aren't quite set up
//   //to integrate to one.
//   return Gradient;
// }
//
// Eigen::MatrixXd PolynomialChaosExpansion::ComputeVarianceHessian() const
// {
//   VectorXd normalVec = GetNormalizationVec();
//
//   //if there's only the one constant term, the PCE has variance zero
//   if (normalVec.rows() <= 1) {
//     return MatrixXd::Zero(coeffs.cols(), coeffs.cols());
//   }
//
//   normalVec(0) = 0;
//   VectorXd squareNorm = normalVec.array().square();
//
//   VectorXd dimMeasures(inputSizes(0));
//   for (unsigned int i = 0; i < inputSizes(0); ++i) {
//     dimMeasures(i) = polys[i]->GetNormalization(0);
//   }
//
//   return ((squareNorm / dimMeasures.prod())).asDiagonal();
// }

Eigen::VectorXd PolynomialChaosExpansion::Magnitude() const
{
  Eigen::VectorXd normalVec = GetNormalizationVec();

  //take the matrix product between the squared coeffs and squared norms, then sqrt each term
  return (coeffs.array().square().matrix() * normalVec.array().square().matrix()).array().sqrt();
}

// PolynomialChaosExpansion::Ptr PolynomialChaosExpansion::ComputeWeightedSum(
//   std::vector<PolynomialChaosExpansion::Ptr> expansions,
//   VectorXd const& weights,
//   VariableCollection::Ptr variables,
//   unsigned int const outputDim,
//   UnstructuredMultiIndexFamily::Ptr const& polynomials,
//   std::map<unsigned int, std::set<unsigned int> > const& basisMultiIndexCache)
// {
//   PolynomialChaosExpansion::Ptr sumPCE(new PolynomialChaosExpansion(variables, outputDim));
//
//   //we already know exactly which terms in the pce we need, it was passed in
//   sumPCE->terms = polynomials;
//
//   //allocate space for all the coeffs
//   sumPCE->coeffs = MatrixXd::Zero(outputDim, sumPCE->terms->GetNumberOfIndices());
//
//   //find the non-zero weights
//   VectorXu nonZeroWeights = FindNonZeroInVector(weights);
//
//   //loop over every polynomial in the result
//   for (unsigned int i = 0; i < sumPCE->terms->GetNumberOfIndices(); i++) {
//     //allocate space for the values of the individual terms
//     MatrixXd termCoeffs = MatrixXd::Zero(outputDim, weights.rows());
//
//     //save the multi-index we're looking for
//     RowVectorXi multiToFind = sumPCE->terms->IndexToMulti(i);
//
//     //Find the non-zero weighted polys that actually have term i
//     //Start with this sorting because it's usually very small
//     for (unsigned int j = 0; j < nonZeroWeights.rows(); ++j) {
//       //if this one has the current polynomial
//       PolynomialChaosExpansion::Ptr jthEst = expansions.at(nonZeroWeights(j));
//
//       //do this by checking our cache whether this term is included, if the cache was provided
//       std::map<unsigned int, std::set<unsigned int> >::const_iterator it = basisMultiIndexCache.find(i);
//       if (it->second.find(nonZeroWeights(j)) != it->second.end()) {
//         //copy the coeffs into the results
//         unsigned int index = jthEst->terms->MultiToIndex(multiToFind);
//         termCoeffs.col(nonZeroWeights(j)) = jthEst->coeffs.col(index);
//       }
//     }
//
//     //finish by summing the coeffs with the weights and store them in the result
//     sumPCE->coeffs.col(i) = (termCoeffs * weights);
//   }
//
//   return sumPCE;
// }
//
std::shared_ptr<PolynomialChaosExpansion> PolynomialChaosExpansion::ComputeWeightedSum(std::vector<std::shared_ptr<PolynomialChaosExpansion>>    expansions,
                                                                                       Eigen::VectorXd                                    const& weights,
                                                                                       std::shared_ptr<MultiIndexSet>                     const& polynomials,
                                                                                       std::vector<std::vector<unsigned int>>             const& locToGlob)
{
  // Loop over each expansion and add the weighted coefficients
  Eigen::MatrixXd newCoeffs = Eigen::MatrixXd::Zero(expansions.at(0)->outputSizes(0), polynomials->Size());

  for(unsigned int i=0; i<expansions.size(); ++i){
    if(std::abs(weights(i))>10.0*std::numeric_limits<double>::epsilon()){
      for(unsigned int k=0; k<expansions.at(i)->coeffs.cols(); ++k){
        newCoeffs.col(locToGlob.at(i).at(k)) += weights(i)*expansions.at(i)->coeffs.col(k);
      }
    }
  }

  return std::make_shared<PolynomialChaosExpansion>(expansions.at(0)->basisComps, polynomials, newCoeffs);
}

std::shared_ptr<PolynomialChaosExpansion> PolynomialChaosExpansion::ComputeWeightedSum(std::vector<std::shared_ptr<PolynomialChaosExpansion>> expansions,
                                                                                       Eigen::VectorXd                                 const& weights)
{
  assert(weights.size()==expansions.size());

  unsigned int inputDim = expansions.at(0)->inputSizes(0);
  unsigned int outputDim = expansions.at(0)->outputSizes(0);

  // check the output size of each expansion
  for(int i=1; i<expansions.size(); ++i){
    assert(outputDim==expansions.at(0)->outputSizes(0));
    assert(inputDim==expansions.at(0)->inputSizes(0));
  }

  // the collection of all polynomials that will be in the result
  auto allPolynomials = std::make_shared<MultiIndexSet>(inputDim);

  // Maps the local indices in each expansion to the global index of the sum
  std::vector<std::vector<unsigned int>> locToGlob(weights.size());

  // Make a union of all the multiindex sets, and keep track of the coefficient mapping
  for(int j=0; j<expansions.size(); ++j)
  {
    // if the weight is zero, skip this expansion
    if(std::abs(weights(j))>10.0*std::numeric_limits<double>::epsilon()){

      //loop over all the polynomials in this expansion
      unsigned int numTerms = expansions.at(j)->multis->Size();
      locToGlob.at(j).resize(numTerms);
      for (unsigned int i = 0; i < numTerms; ++i)
        locToGlob.at(j).at(i) = allPolynomials->AddActive(expansions.at(j)->multis->IndexToMulti(i));
    }
  }

  return ComputeWeightedSum(expansions,weights,allPolynomials,locToGlob);
}


Eigen::VectorXd PolynomialChaosExpansion::TotalSensitivity(unsigned const targetDim) const
{

  //grab the total normalizations
  Eigen::VectorXd normalVec = GetNormalizationVec();

  // Zero out the terms that do not depend on the targetDim
  for (unsigned int i = 0; i < multis->Size(); ++i)
    normalVec(i) *=  static_cast<double>((multis->IndexToMulti(i)->GetValue(targetDim) != 0));

  Eigen::VectorXd dimMeasures(inputSizes(0));
  for (unsigned int i = 0; i < inputSizes(0); ++i)
    dimMeasures(i) = std::dynamic_pointer_cast<OrthogonalPolynomial>(basisComps.at(i))->Normalization(0);

  //we want the sum of normalized involved coeffs divided by the normalized sum of all of them.
  //Don't forget the constant normalization the variance requires
  return (coeffs.array().square().matrix() *
          normalVec.array().square().matrix()).array() / dimMeasures.prod() / Variance().array();
}

Eigen::MatrixXd PolynomialChaosExpansion::TotalSensitivity() const
{
  Eigen::MatrixXd result(coeffs.rows(), inputSizes(0));

  for (unsigned int i = 0; i < inputSizes(0); ++i)
    result.col(i) = TotalSensitivity(i);

  return result;
}

Eigen::VectorXd PolynomialChaosExpansion::MainSensitivity(unsigned const targetDim) const
{
  //we're going to fill this index with either 1 or 0, depending on whether the term includes the target dimension
  Eigen::ArrayXd contributesToIndex = Eigen::VectorXd::Zero(NumTerms());

  //create a row vector with a one at only the place we're interested in
  Eigen::Array<int, 1, Eigen::Dynamic> deltaDim = Eigen::RowVectorXi::Ones(inputSizes(0));
  deltaDim(targetDim) = 0;

  //each one contributes iff it is the only non-zero one
  for (unsigned int i = 0; i < multis->Size(); ++i) {
    Eigen::RowVectorXi oneIndex = multis->IndexToMulti(i)->GetVector();
    contributesToIndex(i) = static_cast<double>((oneIndex(targetDim) != 0) && ((oneIndex.array() * deltaDim).matrix().sum() == 0));
  }

  //grab the total normalizations
  Eigen::VectorXd normalVec = GetNormalizationVec();

  //filter the normalizations so that the ones that don't involve targetDim are zero
  Eigen::VectorXd filterdVec = normalVec.array() * contributesToIndex;

  Eigen::VectorXd dimMeasures(inputSizes(0));
  for (unsigned int i = 0; i < inputSizes(0); ++i)
    dimMeasures(i) = std::dynamic_pointer_cast<OrthogonalPolynomial>(basisComps.at(i))->Normalization(0);

  //we want the sum of normalized involved coeffs divided by the normalized sum of all of them.
  //Don't forget the constant normalization the variance requires
  return (coeffs.array().square().matrix() *
          filterdVec.array().square().matrix()).array() / dimMeasures.prod() / Variance().array();
}

Eigen::MatrixXd PolynomialChaosExpansion::MainSensitivity() const
{
  Eigen::MatrixXd result(coeffs.rows(), inputSizes(0));

  for (unsigned int i = 0; i < inputSizes(0); ++i)
    result.col(i) = MainSensitivity(i);

  return result;
}
