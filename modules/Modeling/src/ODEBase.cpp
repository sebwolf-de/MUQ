#include "MUQ/Modeling/ODEBase.h"

using namespace muq::Modeling;

ODEBase::ODEBase(std::shared_ptr<WorkPiece> rhs) : WorkPiece(), rhs(rhs) {
  // we must know the number of inputs for the rhs and it must have at least one (the state)
  assert(rhs->numInputs>0);

  // set the input and output types
  SetInputOutputTypes();
}

ODEBase::~ODEBase() {}

void ODEBase::SetInputOutputTypes() {
  // the type of the first input (the state) for the rhs and the root
  const std::string& rhsStateType = rhs->InputType(0, false);


  // the first input and output type is the state type --- if the type is known the rhs and the root must agree
  if( rhsStateType.compare("")!=0 ) { // we know the state type for the rhs
    inputTypes[0] = rhsStateType;
    outputTypes[0] = rhsStateType;
  }

  // the next set of input parameters are the parameters for the rhs
  for( auto intype : rhs->InputTypes() ) {
    if( intype.first==0 ) { // we've already set the state type
      continue;
    }

    inputTypes[intype.first] = intype.second;
  }
}

bool ODEBase::CheckFlag(void* flagvalue, std::string const& funcname, unsigned int const opt) const {
  // there are only two options
  assert(opt==0 || opt==1);

  // check if Sundials function returned nullptr pointer - no memory allocated 
  if( opt==0 && flagvalue==nullptr ) { 
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname.c_str());

    // return failure 
    return false; 
  }

  // check if flag<0 
  if( opt==1 ) {
    // get int value
    int *errflag = (int *) flagvalue;
    
    if( *errflag<0 ) { // negative indicates Sundials error
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname.c_str(), *errflag);

      // return failure
      return false; 
    }
  }
  
  // return success
  return true;
}
