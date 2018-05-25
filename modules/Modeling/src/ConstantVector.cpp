#include "MUQ/Modeling/ConstantVector.h"
#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Modeling;


ConstantVector::ConstantVector(Eigen::VectorXd const& valIn) : ModPiece(Eigen::VectorXi(), valIn.size()*Eigen::VectorXi::Ones(1))
{
  outputs.resize(1);
  outputs.at(0) = valIn;
}

void ConstantVector::SetValue(Eigen::VectorXd const& valIn)
{
  if(valIn.size() != outputSizes(0))
    throw muq::WrongSizeError("In ConstantVector::SetValue, new vector has size " + std::to_string(valIn.size()) + ", but expected a size of " + std::to_string(outputSizes(0)) + ".");

  outputs.at(0) = valIn;
}

void ConstantVector::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs){
  return;
}
