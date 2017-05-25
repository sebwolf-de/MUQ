#include "MUQ/Modeling/RootfindingIVP.h"

using namespace muq::Modeling;

RootfindingIVP::RootfindingIVP(std::shared_ptr<WorkPiece> rhs, std::shared_ptr<WorkPiece> root) : WorkPiece(), rhs(rhs), root(root) {
  // we must know the number of inputs for both the rhs and the root and they must have at least one (the state)
  assert(rhs->numInputs>0);
  assert(root->numInputs>0);

  // set the input and output types
  SetInputOutputTypes();
}

RootfindingIVP::~RootfindingIVP() {}

void RootfindingIVP::EvaluateImpl(ref_vector<boost::any> const& inputs) {
}

void RootfindingIVP::SetInputOutputTypes() {
  // the type of the first input (the state) for the rhs and the root
  const std::string& rhsStateType = rhs->InputType(0, false);
  const std::string& rootStateType = root->InputType(0, false);

  // the first input and output type is the state type --- if the type is known the rhs and the root must agree
  assert(rhsStateType.compare("")==0 || rootStateType.compare("")==0 || rootStateType.compare(rhsStateType)==0);
  if( rhsStateType.compare("")!=0 ) { // we know the state type for the rhs
    inputTypes[0] = rhsStateType;
    outputTypes[0] = rhsStateType;
  } else if( rootStateType.compare("")!=0 ) { // we don't know the state type for the rhs, but do for the root
    inputTypes[0] = rootStateType;
    outputTypes[0] = rootStateType;
  }

  // the next set of input parameters are the parameters for the rhs
  for( auto intype : rhs->InputTypes() ) {
    if( intype.first==0 ) { // we've already set the state type
      continue;
    }

    inputTypes[intype.first] = intype.second;
  }

  // the next set of input parameters are the parameters for the root
  for( auto intype : root->InputTypes() ) {
    if( intype.first==0 ) { // we've already set the state type
      continue;
    }

    inputTypes[rhs->numInputs+intype.first-1] = intype.second;
  }
}
