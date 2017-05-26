#include "MUQ/Modeling/RootfindingIVP.h"

// Sundials includes
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>

using namespace muq::Modeling;

RootfindingIVP::RootfindingIVP(std::shared_ptr<WorkPiece> rhs, std::shared_ptr<WorkPiece> root, std::shared_ptr<AnyAlgebra> algebra) : ODEBase(rhs), root(root), algebra(algebra) {
  // we must know the number of inputs for both the rhs and the root and they must have at least one (the state)
  assert(rhs->numInputs>0);
  assert(root->numInputs>0);

  // set the input and output types
  UpdateInputOutputTypes();
}

RootfindingIVP::~RootfindingIVP() {}

void RootfindingIVP::CVODES(ref_vector<boost::any> const& inputs) const {
  // get the dimension of the state
  const unsigned int dim = algebra->VectorDimensionBase(inputs[0]);
  assert(dim>0);
  
  // initialize the state (set the initial conditions)
  N_Vector state = nullptr;
  InitializeState(state, inputs[0], dim);
}

void RootfindingIVP::InitializeState(N_Vector& state, boost::any const& ic, unsigned int const dim) const {
  // initialize the state (set the initial conditions)
  state = N_VNew_Serial(dim);
  assert(CheckFlag((void*)state, "N_VNew_Serial", 0)); // make sure state was properly initialized

  // set the values to the initial conditions
  for( unsigned int i=0; i<dim; ++i ) {
    // NV_Ith_S references the ith component of the vector v
    NV_Ith_S(state, i) = boost::any_cast<double>(algebra->AccessElementBase(i, ic));
  }
}

void RootfindingIVP::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  CVODES(inputs);
}

void RootfindingIVP::UpdateInputOutputTypes() {
  // the type of the first input (the state) for the rhs and the root
  const std::string& rhsStateType = rhs->InputType(0, false);
  const std::string& rootStateType = root->InputType(0, false);

  // the first input and output type is the state type --- if the type is known the rhs and the root must agree
  assert(rhsStateType.compare("")==0 || rootStateType.compare("")==0 || rootStateType.compare(rhsStateType)==0);
  if( rhsStateType.compare("") && rootStateType.compare("")!=0 ) { // we don't know the state type for the rhs, but do for the root
    inputTypes[0] = rootStateType;
    outputTypes[0] = rootStateType;
  }

  // the second set of input parameters are the parameters for the root
  for( auto intype : root->InputTypes() ) {
    if( intype.first==0 ) { // we've already set the state type
      continue;
    }

    inputTypes[rhs->numInputs+intype.first-1] = intype.second;
  }

  // the second output is the time where the root is reached
  outputTypes[1] = typeid(double).name();
}
