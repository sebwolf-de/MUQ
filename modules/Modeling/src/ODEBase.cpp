#include "MUQ/Modeling/ODEBase.h"

#include <cvodes/cvodes.h> // prototypes for CVODE fcts. and consts. 
#include <cvodes/cvodes_spgmr.h> // prototypes & constants for CVSPGMR solver 
#include <cvodes/cvodes_spbcgs.h> // prototypes & constants for CVSPBCG solver 
#include <cvodes/cvodes_sptfqmr.h> // prototypes & constants for SPTFQMR solver 
#include <cvodes/cvodes_dense.h> // prototype for CVDense 

#include <sundials/sundials_dense.h> // definitions DlsMat DENSE_ELEM 
#include <sundials/sundials_types.h> // definition of type 
#include <sundials/sundials_math.h>  // contains the macros ABS, SQR, and EXP 

namespace pt = boost::property_tree;
using namespace muq::Modeling;

ODEBase::ODEBase(std::shared_ptr<WorkPiece> rhs, pt::ptree const& pt, std::shared_ptr<AnyAlgebra> algebra) : WorkPiece(), rhs(rhs), algebra(algebra), multiStep(pt.get<std::string>("ODESolver.MultistepMethod", "BDF")), nonlinSolver(pt.get<std::string>("ODESolver.Solver", "Newton")) {
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

void ODEBase::InitializeState(N_Vector& state, boost::any const& ic, unsigned int const dim) const {
  // initialize the state (set the initial conditions)
  state = N_VNew_Serial(dim);
  assert(CheckFlag((void*)state, "N_VNew_Serial", 0)); // make sure state was properly initialized

  // set the values to the initial conditions
  for( unsigned int i=0; i<dim; ++i ) {
    // NV_Ith_S references the ith component of the vector v
    NV_Ith_S(state, i) = boost::any_cast<double>(algebra->AccessElementBase(i, ic));
  }
}

void ODEBase::CreateSolverMemory(void* cvode_mem) const {
  // determine the multistep method and the nonlinear solver
  assert(multiStep.compare("Adams")==0 || multiStep.compare("BDF")==0);
  const int multMethod = (multiStep.compare("BDF")==0) ? CV_BDF : CV_ADAMS;
  assert(nonlinSolver.compare("Iter")==0 || nonlinSolver.compare("Newton")==0);
  const int solveMethod = (nonlinSolver.compare("Newton")==0) ? CV_NEWTON : CV_FUNCTIONAL;

  // create the memory
  cvode_mem = CVodeCreate(multMethod, solveMethod);
  assert(CheckFlag((void*)cvode_mem, "CVodeCreate", 0));
}
