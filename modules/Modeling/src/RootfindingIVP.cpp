#include "MUQ/Modeling/RootfindingIVP.h"

// Sundials includes
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <cvodes/cvodes_dense.h>     /* prototype for CVDense */

namespace pt = boost::property_tree;
using namespace muq::Modeling;

RootfindingIVP::RootfindingIVP(std::shared_ptr<WorkPiece> rhs, std::shared_ptr<WorkPiece> root, pt::ptree const& pt, std::shared_ptr<AnyAlgebra> algebra) : ODEBase(rhs, pt, algebra), root(root) {
  // we must know the number of inputs for both the rhs and the root and they must have at least one (the state)
  assert(rhs->numInputs>0);
  assert(root->numInputs>0);

  // we must know the number of outputs for the root
  assert(root->numOutputs>0);

  // set the input and output types
  UpdateInputOutputTypes();
}

RootfindingIVP::~RootfindingIVP() {}

void RootfindingIVP::FindRoot(ref_vector<boost::any> const& inputs) {
  // the number of inputs must be at leastthan the number of inputs required by the rhs and the root
  assert(inputs.size()>=rhs->numInputs+root->numInputs-1);

  // create the state vector (have to do a hard copy --- N_Vector is a pointer to the data, the pointer has been declared const, not the data)
  N_Vector state;
  DeepCopy(state, boost::any_cast<const N_Vector&>(inputs[0]));

  // create a data structure to pass around in Sundials
  auto data = std::make_shared<ODEData>(rhs, root, inputs);

  // set solver to null
  void* cvode_mem = nullptr;

  // create the solver memory
  cvode_mem = CVodeCreate(multiStep, solveMethod);
  assert(CheckFlag((void*)cvode_mem, "CVodeCreate", 0));

  // initialize the solver
  CreateSolverMemory(cvode_mem, state, data);

  // tell the solver how evaluate the root
  int flag = CVodeRootInit(cvode_mem, root->numOutputs, EvaluateRoot);
  assert(CheckFlag(&flag, "CVodeRootInit", 1));

  realtype t=0.0;
  flag = CVode(cvode_mem, 10.0, state, &t, CV_NORMAL);

  // set the output
  outputs.push_back(state);
}

void RootfindingIVP::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  FindRoot(inputs);
}

void RootfindingIVP::UpdateInputOutputTypes() {
  // the type of the first input (the state) for the rhs and the root
  assert(rhs->InputType(0, false).compare(root->InputType(0, false))==0);

  // the second set of input parameters are the parameters for the root
  for( auto intype : root->InputTypes() ) {
    if( intype.first==0 ) { // we've already set the state type
      continue;
    }

    inputTypes[rhs->numInputs+intype.first-1] = intype.second;
  }

  // the first output type is the state
  outputTypes[0] = rhs->InputType(0, false);

  // the second output is the time where the root is reached
  outputTypes[1] = typeid(double).name();
}

int RootfindingIVP::EvaluateRoot(realtype t, N_Vector state, realtype *root, void *user_data) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->root);

  // set the state input
  const boost::any& anyref = state;

  // the inputs the root function
  ref_vector<boost::any> rootins(data->inputs.begin()+data->rhs->numInputs, data->inputs.end());
  rootins.insert(rootins.begin(), anyref);

  // evaluate the root
  const std::vector<boost::any>& result = data->root->Evaluate(rootins);
  assert(result.size()==data->root->numOutputs);
  for( unsigned int i=0; i<result.size(); ++i ) {
    root[i] = boost::any_cast<double>(result[i]);
  }

  return 0;
}
