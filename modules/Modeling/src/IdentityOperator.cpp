#include "MUQ/Modeling/IdentityOperator.h"
#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Modeling;


IdentityOperator::IdentityOperator(unsigned int dim) : ModPiece(int(dim)*Eigen::VectorXi::Ones(1),
                                                                int(dim)*Eigen::VectorXi::Ones(1))
{}

void IdentityOperator::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs){
  outputs.resize(1);
  outputs.at(0) = inputs.at(0);
}
