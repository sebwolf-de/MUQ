#include "MUQ/Approximation/Quadrature/SmolyakQuadrature.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"

#include <unordered_map>

using namespace muq::Approximation;
using namespace muq::Utilities;

SmolyakQuadrature::SmolyakQuadrature(std::vector<std::shared_ptr<Quadrature>> scalarRulesIn) : Quadrature(scalarRulesIn.size()),
                                                                                               scalarRules(scalarRulesIn)
{}


void SmolyakQuadrature::Compute(unsigned int order)
{
  Compute(order*Eigen::RowVectorXi::Ones(dim));
}

void SmolyakQuadrature::Compute(Eigen::RowVectorXi const& orders)
{

  // Build the multiindex set.
  std::shared_ptr<MultiIndexSet> multis = BuildMultis(orders);

  Compute(multis);
}


std::shared_ptr<MultiIndexSet> SmolyakQuadrature::BuildMultis(Eigen::RowVectorXi const& orders) const
{
  // Use the minimum order in the vector to form a total order multiindex
  int minOrder = orders.minCoeff();
  assert(minOrder>=0);

  auto multis = MultiIndexFactory::CreateTotalOrder(dim,minOrder);

  // Add other terms to get the right order
  for(int i=0; i<dim; ++i){
    for(int p=minOrder+1; p<=orders(i); ++p)
    {
      auto newMulti = std::make_shared<MultiIndex>(dim);
      newMulti->SetValue(i,p);
      multis += newMulti;
    }
  }

  return multis;
}


void SmolyakQuadrature::Compute(std::shared_ptr<MultiIndexSet> const& multis) {

  // Compute the weights caused by using a tensor product of quadrature rules directly
  Eigen::VectorXd smolyWts = ComputeWeights(multis);

  auto tensorQuad = std::make_shared<FullTensorQuadrature>(scalarRules);

  // A map holding the unique quadrature points (keys) and weights (values)
  std::unordered_map<Eigen::VectorXd, double> quadParts;

  for(int i=0; i<multis->Size(); ++i)
  {
    Eigen::RowVectorXi multiVec = multis->IndexToMulti(i)->GetVector();
    tensorQuad->Compute(multiVec);

    auto& tensorPts = tensorQuad->Points();
    auto& tensorWts = tensorQuad->Weights();

    // Add all the tensor product points to the unordered_map of Smolyak points
    for(int ptInd = 0; ptInd<tensorPts.cols(); ++ptInd){

      auto iter = quadParts.find(tensorPts.col(ptInd));
      if(iter!=quadParts.end()){
        iter->second += smolyWts(i)*tensorWts(ptInd);
      }else{
        quadParts[tensorPts.col(ptInd)] = smolyWts(i)*tensorWts(ptInd);
      }
    }

  }

  // Copy the unordered map, which has unique keys, into the Eigen matrix
  pts.resize(dim,quadParts.size());
  wts.resize(quadParts.size());
  unsigned int ind = 0;
  for(auto& part : quadParts){
    pts.col(ind) = part.first;
    wts(ind) = part.second;
    ++ind;
  }
}


Eigen::VectorXd SmolyakQuadrature::ComputeWeights(std::shared_ptr<MultiIndexSet> const& multis) const
{
  Eigen::VectorXd multiWeights = Eigen::VectorXd::Zero(multis->Size());

  for(unsigned int i = 0; i<multis->Size(); ++i) {
    auto& k = multis->IndexToMulti(i);

    multiWeights(i) += 1.0;

    for(unsigned int d = 0; d<multis->GetMultiLength(); ++d)
    {
      auto newMulti = std::make_shared<MultiIndex>(k);
      unsigned int oldVal = k->GetValue(d);
      if(k->GetValue(d)>0){
        newMulti->SetValue(d,oldVal-1);
        multiWeights(multis->MultiToIndex(newMulti)) -= 1.0;
      }
    }
  }

  return multiWeights;
}
