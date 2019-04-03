#include "MUQ/Approximation/PolynomialChaos/PCEFactory.h"
#include "MUQ/Approximation/Quadrature/FullTensorQuadrature.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

using namespace muq::Approximation;
using namespace muq::Modeling;
using namespace muq::Utilities;

PCEFactory::PCEFactory(std::vector<std::shared_ptr<Quadrature>>         const& quadTypes,
                       std::shared_ptr<muq::Utilities::MultiIndex>      const& quadOrders,
                       std::vector<std::shared_ptr<IndexedScalarBasis>> const& polyTypesIn) : polyTypes(polyTypesIn)
{
  unsigned int dim = quadTypes.size();
  assert(dim==quadOrders->GetLength());
  assert(dim==polyTypesIn.size());

  /* Create a multiindex set containing polynomials that can be integrated exactly
     with a tensor product quadrature rule defined by quadOrders and quadTypesIn
  */
  Eigen::RowVectorXi orders(dim);
  for(unsigned int i=0; i<dim; ++i)
    orders(i) = std::floor( quadTypes.at(i)->Exactness(quadOrders->GetValue(i)) / 2.0);

  polyMultis = MultiIndexFactory::CreateFullTensor(orders);

  // Build the tensor product quadrature rule
  FullTensorQuadrature tensQuad(quadTypes, quadOrders->GetVector());
  quadWts = tensQuad.Weights();
  quadPts.resize( tensQuad.Points().cols());
  for(unsigned int i=0; i<quadPts.size(); ++i)
    quadPts.at(i) = tensQuad.Points().col(i);
}


PCEFactory::PCEFactory(std::vector<std::shared_ptr<Quadrature>>          const& quadTypes,
                       std::shared_ptr<muq::Utilities::MultiIndex>       const& quadOrders,
                       std::vector<std::shared_ptr<IndexedScalarBasis>>  const& polyTypesIn,
                       std::shared_ptr<MultiIndexSet>                    const& polyMultisIn) : polyTypes(polyTypesIn),
                                                                                                polyMultis(polyMultisIn)
{
  unsigned int dim = quadTypes.size();
  assert(dim==quadOrders->GetLength());
  assert(dim==polyTypesIn.size());
  assert(dim==polyMultisIn->GetMultiLength());

  // Build the tensor product quadrature rule
  FullTensorQuadrature tensQuad(quadTypes, quadOrders->GetVector());
  quadWts = tensQuad.Weights();
  quadPts.resize( tensQuad.Points().cols());
  for(unsigned int i=0; i<quadPts.size(); ++i)
    quadPts.at(i) = tensQuad.Points().col(i);
}

std::shared_ptr<PolynomialChaosExpansion> PCEFactory::Compute(std::shared_ptr<muq::Modeling::ModPiece> const& model)
{
  // Make sure the model has a single input and single output
  assert(model->outputSizes.size()==1);
  assert(model->inputSizes.size()==1);

  // Evaluate the model at the quadrature points
  std::vector<Eigen::VectorXd> quadEvals(quadPts.size());
  for(unsigned int k=0; k<quadPts.size(); ++k)
    quadEvals.at(k) = model->Evaluate(quadEvals.at(k)).at(0);

  return Compute(quadEvals);
}

std::shared_ptr<PolynomialChaosExpansion> PCEFactory::Compute(std::vector<Eigen::VectorXd> const& quadEvals)
{
  assert(quadEvals.size()==quadPts.size());
  unsigned int outputDim = quadEvals.at(0).size();

  auto pce = std::make_shared<PolynomialChaosExpansion>(polyTypes, polyMultis, outputDim);

  // Compute all of the coefficients
  Eigen::MatrixXd coeffs = Eigen::MatrixXd::Zero(quadEvals.at(0).size(), polyMultis->Size());

  // Estimate the projection c_{ij} = <\psi_j, f_i> = \sum_{k=1}^N w_k \psi_j(x_k) f_i(x_k)
  for(unsigned int k=0; k<quadPts.size(); ++k)
    coeffs += quadWts(k) * quadEvals.at(k) * pce->BuildVandermonde(quadPts.at(k));

  Eigen::VectorXd allNorms = pce->GetNormalizationVec().array().square().matrix();
  coeffs = (coeffs.array().rowwise() / allNorms.transpose().array()).eval();

  pce->SetCoeffs(coeffs);

  return pce;
}
