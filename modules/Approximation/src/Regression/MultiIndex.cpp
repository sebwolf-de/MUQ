#include "MUQ/Approximation/Regression/MultiIndex.h"

using namespace muq::Modeling;
using namespace muq::Approximation;

MultiIndex::MultiIndex(unsigned int const n, unsigned int const maxOrder, unsigned int const minOrder) :
  WorkPiece(std::vector<std::string>({typeid(unsigned int).name()}), std::vector<std::string>({typeid(Eigen::VectorXi).name()})) {
  // add all of the mutli indices that satify the order conditions
  alphas.clear();
  AddMulti(Eigen::VectorXi::Zero(n), minOrder, maxOrder);
}

MultiIndex::~MultiIndex() {}

void MultiIndex::AddMulti(Eigen::VectorXi const alpha, unsigned int const minOrder, unsigned int const maxOrder, unsigned int const dim) {
  // compute the order of this multi index
  unsigned int const order = alpha.sum();

  // if the order is too large, return
  if( order>maxOrder ) { return; }

  // we know the order is not to large---if the order is also not too small, save it
  if( minOrder<=order ) { alphas.push_back(alpha); }

  // loop through the dimension
  for( unsigned int i=dim; i<alpha.size(); ++i ) {
    // create a new multi-index adding one to the i^{th} dimension
    Eigen::VectorXi alpha_new = alpha;
    alpha_new(i) += 1;

    // try adding the new alpha
    AddMulti(alpha_new, minOrder, maxOrder, i);
  }
}

unsigned int MultiIndex::Size() const {
  return alphas.size();
}

void MultiIndex::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // determine which multi-index we need
  const unsigned int i = boost::any_cast<unsigned int const>(inputs[0]);
  assert(i<alphas.size());
  
  // retrive the ith multi-index
  outputs.resize(1);
  outputs[0] = alphas[i];
}
