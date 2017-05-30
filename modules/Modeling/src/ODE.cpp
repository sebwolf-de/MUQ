#include "MUQ/Modeling/ODE.h"

#include <cvodes/cvodes.h> // prototypes for CVODE fcts. and consts. 
#include <cvodes/cvodes_spgmr.h> // prototypes & constants for CVSPGMR solver 
#include <cvodes/cvodes_spbcgs.h> // prototypes & constants for CVSPBCG solver 
#include <cvodes/cvodes_sptfqmr.h> // prototypes & constants for SPTFQMR solver 
#include <cvodes/cvodes_dense.h> // prototype for CVDense 

#include <sundials/sundials_dense.h> /* use generic DENSE solver in preconditioning */
#include <sundials/sundials_types.h> // definition of type 
#include <sundials/sundials_math.h>  // contains the macros ABS, SQR, and EXP 

namespace pt = boost::property_tree;
using namespace muq::Modeling;

ODE::ODE(std::shared_ptr<WorkPiece> rhs, pt::ptree const& pt, std::shared_ptr<AnyAlgebra> algebra) : ODEBase(rhs, pt, algebra) {}

ODE::~ODE() {}

void ODE::Integrate(ref_vector<boost::any> const& inputs) {
  // create the state vector (have to do a hard copy --- N_Vector is a pointer to the data, the pointer has been declared const, not the data)
  const N_Vector& ic = boost::any_cast<N_Vector>(inputs[0]);
  N_Vector state = N_VNew_Serial(NV_LENGTH_S(ic));
  for( unsigned int i=0; i<NV_LENGTH_S(ic); ++i ) {
    NV_Ith_S(state, i) = NV_Ith_S(ic, i);
  }

  // create a data structure to pass around in Sundials
  auto data = std::make_shared<ODEData>(rhs, inputs);

  // set solver to null
  void* cvode_mem = nullptr;

  // create the solver memory
  cvode_mem = CVodeCreate(multiStep, solveMethod);
  assert(CheckFlag((void*)cvode_mem, "CVodeCreate", 0));

  // initialize the solver
  CreateSolverMemory(cvode_mem, state, data);
  
  // each element corresponds to a vector of desired times, first: the current index of that vector, second: the size of that vector
  std::vector<std::pair<unsigned int, unsigned int> > timeIndices(inputs.size()-rhs->numInputs);
  const boost::any none = boost::none;
  outputs.resize(timeIndices.size(), none);
  for( unsigned int i=0; i<timeIndices.size(); ++i ) { // loop through the desired outputs
    timeIndices[i].first = 0; // start at the first index
    timeIndices[i].second = algebra->VectorDimensionBase(inputs[rhs->numInputs+i]); // the size of the vector

    if( timeIndices[i].second>1 ) {
      outputs[i] = std::vector<N_Vector>();
      std::vector<N_Vector>& outref = boost::any_cast<std::vector<N_Vector>&>(outputs[i]);
      outref.reserve(timeIndices[i].second);
    }
  }

  // first: the next time to integrate to, second: the output index
  std::pair<double, int> nextTime;

  double tcurrent = 0.0;
  while( NextTime(nextTime, timeIndices, ref_vector<boost::any>(inputs.begin()+rhs->numInputs, inputs.end())) ) {
    if( fabs(nextTime.first-tcurrent)>1.0e-14 ) { // we have to move forward in time
      int flag = CVode(cvode_mem, nextTime.first, state, &tcurrent, CV_NORMAL);
      assert(CheckFlag(&flag, "CVode", 1));

      if( timeIndices[nextTime.second].second>1 ) {
	std::vector<N_Vector>& outref = boost::any_cast<std::vector<N_Vector>&>(outputs[nextTime.second]);

	outref.push_back(N_VNew_Serial(NV_LENGTH_S(state)));
	for( unsigned int i=0; i<NV_LENGTH_S(state); ++i ) {
	  NV_Ith_S(outref[outref.size()-1], i) = NV_Ith_S(state, i);
	}
      } else {
	outputs[nextTime.second] = N_VNew_Serial(NV_LENGTH_S(state));
	N_Vector& outref = boost::any_cast<N_Vector&>(outputs[nextTime.second]);

	for( unsigned int i=0; i<NV_LENGTH_S(state); ++i ) {
	  NV_Ith_S(outref, i) = NV_Ith_S(state, i);
	}
      }
    } else {
      if( timeIndices[nextTime.second].second>1 ) {
	std::vector<N_Vector>& outref = boost::any_cast<std::vector<N_Vector>&>(outputs[nextTime.second]);

	outref.push_back(N_VNew_Serial(NV_LENGTH_S(state)));
	for( unsigned int i=0; i<NV_LENGTH_S(state); ++i ) {
	  NV_Ith_S(outref[outref.size()-1], i) = NV_Ith_S(state, i);
	}
      } else {
	outputs[nextTime.second] = N_VNew_Serial(NV_LENGTH_S(state));
	N_Vector& outref = boost::any_cast<N_Vector&>(outputs[nextTime.second]);

	for( unsigned int i=0; i<NV_LENGTH_S(state); ++i ) {
	  NV_Ith_S(outref, i) = NV_Ith_S(state, i);
	}
      }
    }
  }


  //realtype t=0.0;
  //flag = CVode(cvode_mem, 1.0, state, &t, CV_NORMAL);
  //PrintOutput(t, NV_Ith_S(state, 0), NV_Ith_S(state, 1));

  /*// the number of inputs must be greater than the number of inputs required by the rhs
  assert(inputs.size()>rhs->numInputs);
  
  // get the dimension of the state
  const unsigned int dim = algebra->VectorDimensionBase(inputs[0]);
  assert(dim>0);

  // the evolving state
  N_Vector state = nullptr;
  state = N_VNew_Serial(dim);
  state = boost::any_cast<N_Vector>(inputs[0]);
  assert(CheckFlag((void *)state,  "N_VNew_Serial", 0));
  std::cout << "INITIAL STATE: " << NV_Ith_S(state, 0) << " " << NV_Ith_S(state, 1) << std::endl << std::endl;
   
  // create a data structure to pass around in Sundials
  auto data = std::make_shared<ODEData>(rhs, ref_vector<boost::any>(inputs.begin(), inputs.begin()+rhs->numInputs));

  // crete the solver memory
  void* cvode_mem = nullptr;

  // determine the multistep method and the nonlinear solver
  assert(multiStep.compare("Adams")==0 || multiStep.compare("BDF")==0);
  const int multMethod = (multiStep.compare("BDF")==0) ? CV_BDF : CV_ADAMS;
  assert(nonlinSolver.compare("Iter")==0 || nonlinSolver.compare("Newton")==0);
  const int solveMethod = (nonlinSolver.compare("Newton")==0) ? CV_NEWTON : CV_FUNCTIONAL;

  // create the memory
  cvode_mem = CVodeCreate(multMethod, solveMethod);
  assert(CheckFlag((void*)cvode_mem, "CVodeCreate", 0));

  // set the pointer to user-defined data
  assert(data);
  int flag = CVodeSetUserData(cvode_mem, data.get());
  assert(CheckFlag(&flag, "CVodeSetUserData", 1));

  // set the pointer to user-defined data 
  flag = CVodeSetErrHandlerFn(cvode_mem, ErrorHandler, data.get());
  assert(CheckFlag(&flag, "CVodeSetErrHandlerFun", 1));

  // tell the solver how to evaluate the rhs
  flag = CVodeInit(cvode_mem, EvaluateRHS, 0.0, state);
  assert(CheckFlag(&flag, "CVodeInit", 1));
  
  // initialize the solver
  //InitializeSolver(cvode_mem, state);

  // call CVodeSStolerances to specify the relative and absolute tolerances
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  assert(CheckFlag(&flag, "CVodeSStolerances", 1));

  // set the maximum time step size
  flag = CVodeSetMaxStep(cvode_mem, maxStepSize);
  assert(CheckFlag(&flag, "CVodeSetMaxStep", 1));

  if( linSolver.compare("Dense")==0 ) {
    // Call CVDense to specify the CVDENSE dense linear solver 
    flag = CVDense(cvode_mem, NV_LENGTH_S(state));
    assert(CheckFlag(&flag, "CVDense", 1));

    // Set the Jacobian routine to Jac (user-supplied) 
    flag = CVDlsSetDenseJacFn(cvode_mem, RHSJacobian);
    assert(CheckFlag(&flag, "CVDlsSetDenseJacFn", 1));
  } else {
    if( linSolver.compare("SPGMR")==0 ) {
      flag = CVSpgmr(cvode_mem, 0, 0);
      assert(CheckFlag(&flag, "CVSpgmr", 1));
    } else if( linSolver.compare("SPBCG")==0 ) {
      flag = CVSpbcg(cvode_mem, 0, 0);
      assert(CheckFlag(&flag, "CVSpbcg", 1));
    } else if( linSolver.compare("SPTFQMR")==0 ) {
      flag = CVSptfqmr(cvode_mem, 0, 0);
      assert(CheckFlag(&flag, "CVSptfqmr", 1));
    } else {
      std::cerr << "\nInvalid CVODES linear solver type.  Options are Dense, SPGMR, SPBCG, or SPTFQMR\n\n";
      assert(false);
    }

    // set the Jacobian-times-vector function 
    flag = CVSpilsSetJacTimesVecFn(cvode_mem, RHSJacobianAction);
    assert(CheckFlag(&flag, "CVSpilsSetJacTimesVecFn", 1));
    }*/

  /*//CreateSolverMemory(cvode_mem, state, data);

  // each element corresponds to a vector of desired times, first: the current index of that vector, second: the size of that vector
  std::vector<std::pair<unsigned int, unsigned int> > timeIndices(inputs.size()-rhs->numInputs);
  for( unsigned int i=0; i<timeIndices.size(); ++i ) { // loop through the desired outputs
    timeIndices[i].first = 0; // start at the first index
    timeIndices[i].second = algebra->VectorDimensionBase(inputs[rhs->numInputs+i]); // the size of the vector
  }

  // first: the next time to integrate to, second: the output index
  std::pair<double, int> nextTime;*/

  // the current time
  //double tcurrent = 0.0;

  // while we should keep integrating
  //flag = CVode(cvode_mem, 0.5, state, &tcurrent, CV_NORMAL);
  //std::cout << flag << std::endl;
  //std::cout << "state: " << NV_Ith_S(state, 0) << " " << NV_Ith_S(state, 1) << std::endl;
  /*while( NextTime(nextTime, timeIndices, ref_vector<boost::any>(inputs.begin()+rhs->numInputs, inputs.end())) ) {
    if( fabs(nextTime.first-tcurrent)>1.0e-14 ) { // we have to move forward in time
      std::cout << "next time: " << nextTime.first << std::endl;
      flag = CVode(cvode_mem, nextTime.first, state, &tcurrent, CV_NORMAL);
      std::cout << flag << std::endl;
      std::cout << "state: " << NV_Ith_S(state, 0) << " " << NV_Ith_S(state, 1) << std::endl;

      std::cout << std::endl;
    } else {
    }
    }*/
}

bool ODE::NextTime(std::pair<double, int>& nextTime, std::vector<std::pair<unsigned int, unsigned int> >& timeIndices, ref_vector<boost::any> const& outputTimes) const {
  // the indices must each correspond to a vector of output times
  assert(timeIndices.size()==outputTimes.size());

  // the next time is inifity and the correspond time vector is invalid
  nextTime.first = std::numeric_limits<double>::infinity();
  nextTime.second = -1;

  // loop through each time vector
  for( unsigned int i=0; i<timeIndices.size(); ++i ) {
    // if we have already computed the state at each time
    if( timeIndices[i].first==timeIndices[i].second ) { continue; }

    // the next time at that vector
    const double t = boost::any_cast<double>(algebra->AccessElementBase(timeIndices[i].first, outputTimes[i]));

    if( t<nextTime.first ) { // if it is the smallest so far ...
      // ... it is the next time and save the corresponding time vector
      nextTime.first = t;
      nextTime.second = i;
    }
  }

  // make sure the index is valid
  if( nextTime.second<0 ) {
    // we are done!
    return false;
  }

  // increment the current index of the relavant time vector
  ++timeIndices[nextTime.second].first;

  // continue integrating
  return true;
}

void ODE::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  Integrate(inputs);
}
