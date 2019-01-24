#include "MUQ/Approximation/Quadrature/FullTensorQuadrature.h"

#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"
#include <iostream>

using namespace muq::Approximation;
using namespace muq::Utilities;

FullTensorQuadrature::FullTensorQuadrature(unsigned int                       dim,
                                           std::shared_ptr<Quadrature> const& rule,
                                           unsigned int                       order) : FullTensorQuadrature(std::vector<std::shared_ptr<Quadrature>>(dim,rule),
                                                                                                            order*Eigen::RowVectorXi::Ones(dim)){};

FullTensorQuadrature::FullTensorQuadrature(std::vector<std::shared_ptr<Quadrature>> const& rulesIn,
                                           Eigen::RowVectorXi                       const& orders) : Quadrature(rulesIn.size()),
                                                                                                       rules(rulesIn)
{
  for(int i=0; i<rules.size(); ++i)
    assert(rules.at(i)->Dim()==1);

  assert(rulesIn.size()==orders.size());

  Compute(orders);
}

void FullTensorQuadrature::Compute(unsigned int order)
{
  Compute(order*Eigen::RowVectorXi::Ones(dim));
}

void FullTensorQuadrature::Compute(Eigen::RowVectorXi const& orders) {

  assert(orders.size()==dim);

  // Compute each 1d rule
  for(int i=0; i<dim; ++i)
    rules.at(i)->Compute(orders(i));

  // Extract the number of points in each rule
  Eigen::RowVectorXi numPts(dim);
  for(int i = 0; i<dim; ++i)
    numPts(i) = rules.at(i)->Weights().size()-1;

  // First, compute all of the multiindices
  auto multis = MultiIndexFactory::CreateFullTensor(numPts);

  // Now compute all the terms
  pts.resize(dim,multis->Size());
  wts = Eigen::VectorXd::Ones(multis->Size());

  for(int i=0; i<multis->Size(); ++i){
    auto& multi = multis->IndexToMulti(i);

    for(int d=0; d<dim; ++d){
      wts(i) *= rules.at(d)->Weights()(multi->GetValue(d));
      pts(d,i) = rules.at(d)->Points()(0,multi->GetValue(d));
    }
  }

}
