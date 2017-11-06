#include "MUQ/Modeling/ODEBase.h"

#include <cvodes/cvodes.h> // prototypes for CVODE fcts. and consts. 
#include <cvodes/cvodes_spgmr.h> // prototypes & constants for CVSPGMR solver 
#include <cvodes/cvodes_spbcgs.h> // prototypes & constants for CVSPBCG solver 
#include <cvodes/cvodes_sptfqmr.h> // prototypes & constants for SPTFQMR solver 
#include <cvodes/cvodes_dense.h> // prototype for CVDense 

#include <sundials/sundials_types.h> // definition of type 
#include <sundials/sundials_math.h>  // contains the macros ABS, SQR, and EXP 

namespace pt = boost::property_tree;
using namespace muq::Modeling;

ODEBase::ODEBase(std::shared_ptr<WorkPiece> rhs, pt::ptree const& pt, std::shared_ptr<AnyAlgebra> algebra) : WorkPiece(), rhs(rhs), algebra(algebra), linSolver(pt.get<std::string>("ODESolver.LinearSolver", "Dense")), reltol(pt.get<double>("ODESolver.RelativeTolerance", 1.0e-8)), abstol(pt.get<double>("ODESolver.AbsoluteTolerance", 1.0e-8)), maxStepSize(pt.get<double>("ODESolver.MaxStepSize", 1.0)), autonomous(pt.get<bool>("ODESolver.Autonomous", true)) {
  // do not clear the outputs --- they have to be destroyed properly
  clearOutputs = false;
  clearDerivatives = false;
  
  // we must know the number of inputs for the rhs and it must have at least one (the state)
  assert(rhs->numInputs>0);

  // set the input and output types
  SetInputOutputTypes();

  // determine the multistep method and the nonlinear solver
  const std::string& multiStepMethod = pt.get<std::string>("ODESolver.MultistepMethod", "BDF");
  assert(multiStepMethod.compare("Adams")==0 || multiStepMethod.compare("BDF")==0); // Adams or BDF
  multiStep = (multiStepMethod.compare("BDF")==0) ? CV_BDF : CV_ADAMS;
  const std::string& nonlinearSolver = pt.get<std::string>("ODESolver.Solver", "Newton"); 
  assert(nonlinearSolver.compare("Iter")==0 || nonlinearSolver.compare("Newton")==0); // Iter or Newton
  solveMethod = (nonlinearSolver.compare("Newton")==0) ? CV_NEWTON : CV_FUNCTIONAL;

  // determine the method of thelinear solver
  std::string linearSolver = pt.get<std::string>("ODESolver.LinearSolver", "Dense");
  if( linearSolver.compare("Dense")==0 ) {
    slvr = LinearSolver::Dense;
  } else if( linSolver.compare("SPGMR")==0 ) {
    slvr = LinearSolver::SPGMR;
  } else if( linSolver.compare("SPBCG")==0 ) {
    slvr = LinearSolver::SPBCG;
  } else if( linSolver.compare("SPTFQMR")==0 ) {
    slvr = LinearSolver::SPTFQMR;
  } else {
    std::cerr << "\nInvalid CVODES linear solver type.  Options are Dense, SPGMR, SPBCG, or SPTFQMR\n\n";
    assert(false);
  }
}

ODEBase::~ODEBase() {}

void ODEBase::SetInputOutputTypes() {
  // the name of an N_Vector type
  const std::string stateType = typeid(N_Vector).name();

  if( autonomous ) {
    // the type of the first input (the state) for the rhs
    assert(stateType.compare(rhs->InputType(0, false))==0);
  } else { // non-autonomous ...
    // the time is the first input
    const std::string doubleType = typeid(double).name();
    assert(doubleType.compare(rhs->InputType(0, false))==0);
    // the state is the second input
    assert(stateType.compare(rhs->InputType(1, false))==0);
  }

  // the first input and output type is the state type --- if the type is known the rhs and the root must agree
  inputTypes[0] = stateType;

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
    NV_Ith_S(state, i) = boost::any_cast<double>(algebra->AccessElement(ic, i));
  }
}

void ODEBase::ErrorHandler(int error_code, const char *module, const char *function, char *msg, void *user_data) {}

void ODEBase::CreateSolverMemory(void* cvode_mem, N_Vector const& state, std::shared_ptr<ODEData> data) const {
  // a flag used for error checking
  int flag;

  // set the pointer to user-defined data
  assert(data);
  flag = CVodeSetUserData(cvode_mem, data.get());
  assert(CheckFlag(&flag, "CVodeSetUserData", 1));

  // tell the solver how to deal with errors
  flag = CVodeSetErrHandlerFn(cvode_mem, ErrorHandler, data.get());
  assert(CheckFlag(&flag, "CVodeSetErrHandlerFun", 1));

  // tell the solver how to evaluate the rhs
  flag = CVodeInit(cvode_mem, EvaluateRHS, 0.0, state);
  assert(CheckFlag(&flag, "CVodeInit", 1));

  // specify the relative and absolute tolerances
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  assert(CheckFlag(&flag, "CVodeSStolerances", 1));

  // set the maximum time step size
  flag = CVodeSetMaxStep(cvode_mem, maxStepSize);
  assert(CheckFlag(&flag, "CVodeSetMaxStep", 1));

  // determine which linear solver to use
  switch( slvr ) {
  case LinearSolver::Dense: { // dense linear solver
    // specify the dense linear solver 
    flag = CVDense(cvode_mem, NV_LENGTH_S(state));
    assert(CheckFlag(&flag, "CVDense", 1));
    
    // set the Jacobian routine to Jac (user-supplied) 
    flag = CVDlsSetDenseJacFn(cvode_mem, RHSJacobian);
    assert(CheckFlag(&flag, "CVDlsSetDenseJacFn", 1));
    break;
  }
  default: { // sparse linear solver 
    if( slvr==LinearSolver::SPGMR ) {
      flag = CVSpgmr(cvode_mem, 0, 0);
      assert(CheckFlag(&flag, "CVSpgmr", 1));
    } else if( slvr==LinearSolver::SPBCG ) {
      flag = CVSpbcg(cvode_mem, 0, 0);
      assert(CheckFlag(&flag, "CVSpbcg", 1));
    } else if( slvr==LinearSolver::SPTFQMR ) {
      flag = CVSptfqmr(cvode_mem, 0, 0);
      assert(CheckFlag(&flag, "CVSptfqmr", 1));
    } else {
      std::cerr << "\nInvalid CVODES linear solver type.  Options are Dense, SPGMR, SPBCG, or SPTFQMR\n\n";
      assert(false);
    }
    
    // set the Jacobian-times-vector function 
    flag = CVSpilsSetJacTimesVecFn(cvode_mem, RHSJacobianAction);
    assert(CheckFlag(&flag, "CVSpilsSetJacTimesVecFn", 1));
  }
  }
}

int ODEBase::EvaluateRHS(realtype time, N_Vector state, N_Vector statedot, void *user_data) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
  const boost::any& anyref = state;
  if( data->autonomous ) {
    data->inputs[0] = anyref; 
  } else {
    const boost::any t = time;
    data->inputs[0] = t; 
    data->inputs[1] = anyref; 
  }

  // evaluate the rhs
  const std::vector<boost::any>& result = data->rhs->Evaluate(ref_vector<boost::any>(data->inputs.begin(), data->inputs.begin()+data->rhs->numInputs));

  // need to do a deep copy to avoid memory leaks
  const N_Vector& vec = boost::any_cast<const N_Vector&>(result[0]);
  assert(NV_LENGTH_S(statedot)==NV_LENGTH_S(vec));
  for( unsigned int i=0; i<NV_LENGTH_S(statedot); ++i ) {
    NV_Ith_S(statedot, i) = NV_Ith_S(vec, i);
  }
  
  return 0;
}

int ODEBase::RHSJacobianAction(N_Vector v, N_Vector Jv, realtype time, N_Vector state, N_Vector rhs, void *user_data, N_Vector tmp) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
  const boost::any anyref = state;
  if( data->autonomous ) {
    data->inputs[0] = anyref; 
  } else {
    const boost::any& t = time;
    data->inputs[0] = t; 
    data->inputs[1] = anyref; 
  }

  // compute the jacobain wrt the state
  const boost::any& jacobianAction = data->rhs->JacobianAction(data->autonomous? 0 : 1, 0, v, ref_vector<boost::any>(data->inputs.begin(), data->inputs.begin()+data->rhs->numInputs));
  N_Vector const& jacAct = boost::any_cast<const N_Vector&>(jacobianAction);
  for( unsigned int i=0; i<NV_LENGTH_S(jacAct); ++i ) {
    NV_Ith_S(Jv, i) = NV_Ith_S(jacAct, i);
  }

  return 0;
}

int ODEBase::RHSJacobian(long int N, realtype time, N_Vector state, N_Vector rhs, DlsMat jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);
  
  // set the state input
  const boost::any& anyref = state;
  if( data->autonomous ) {
    data->inputs[0] = anyref; 
  } else {
    const boost::any t = time;
    data->inputs[0] = t; 
    data->inputs[1] = anyref; 
  }

  // evaluate the jacobian
  const boost::any& jcbn = data->rhs->Jacobian(data->autonomous? 0 : 1, 0, ref_vector<boost::any>(data->inputs.begin(), data->inputs.begin()+data->rhs->numInputs));
  const DlsMat& jcbnref = boost::any_cast<const DlsMat&>(jcbn);
  //DenseCopy(boost::any_cast<const DlsMat&>(jcbn), jac);
  for( unsigned int i=0; i<jcbnref->M; ++i ) {
    for( unsigned int j=0; j<jcbnref->N; ++j ) {
      DENSE_ELEM(jac, i, j) = DENSE_ELEM(jcbnref, i, j);
    }
  }

  return 0;
}

void ODEBase::DeepCopy(N_Vector& copy, N_Vector const& orig) const {
  // initialize the copy
  copy = N_VNew_Serial(NV_LENGTH_S(orig));

  // copy each value
  for( unsigned int i=0; i<NV_LENGTH_S(orig); ++i ) {
    NV_Ith_S(copy, i) = NV_Ith_S(orig, i);
  }
}

void ODEBase::InitializeDerivative(unsigned int const ntimes, unsigned int const stateSize, unsigned int const paramSize, DerivativeMode const& mode) {
  switch( mode ) {
  case DerivativeMode::Jac: { // initialize the jacobian
    if( ntimes==1 ) { // one output --- matrix
      jacobian = NewDenseMat(stateSize, paramSize);
    } else { // multiple outputs --- vector of matrices
      jacobian = std::vector<DlsMat>(ntimes);
    }
    
    return;
  }
  case DerivativeMode::JacAction: { // initialize the action of the jacobian
    if( ntimes==1 ) { // one output --- vector
      jacobianAction = N_VNew_Serial(stateSize);
    } else { // multiple outputs --- vector of vectors
      jacobianAction = std::vector<N_Vector>(ntimes);
    }

    return;
  }
  case DerivativeMode::JacTransAction: { // initialize the action of the jacobian transpose
    if( ntimes==1 ) { // one output --- vector
      jacobianTransposeAction = N_VNew_Serial(paramSize);
    } else { // multiple outputs --- vector of vectors
      jacobianTransposeAction = std::vector<N_Vector>(ntimes);
    }
    
    return;
  }
  default:
    assert(false);
  }
}

void ODEBase::SetUpSensitivity(void *cvode_mem, unsigned int const paramSize, N_Vector *sensState) const {
  // initialze the forward sensitivity solver
  int flag = CVodeSensInit(cvode_mem, paramSize, CV_SIMULTANEOUS, ForwardSensitivityRHS, sensState);
  assert(CheckFlag(&flag, "CVodeSensInit1", 1));

  // set sensitivity tolerances
  Eigen::VectorXd absTolVec = abstol * Eigen::VectorXd::Ones(paramSize);
  flag = CVodeSensSStolerances(cvode_mem, reltol, absTolVec.data());
  assert(CheckFlag(&flag, "CVodeSensSStolerances", 1));

  // error control strategy should test the sensitivity variables
  flag = CVodeSetSensErrCon(cvode_mem, true);
  assert(CheckFlag(&flag, "CVodeSetSensErrCon", 1));
}

int ODEBase::ForwardSensitivityRHS(int Ns, realtype time, N_Vector y, N_Vector ydot, N_Vector *ys, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
  const boost::any& anyref = y;
  data->inputs[0] = anyref;

  // the derivative of the rhs wrt the state
  const boost::any& dfdy_any = data->rhs->Jacobian(0, 0, ref_vector<boost::any>(data->inputs.begin(), data->inputs.begin()+data->rhs->numInputs));
  const DlsMat& dfdy = boost::any_cast<const DlsMat&>(dfdy_any);

  for( unsigned int i=0; i<Ns; ++i ) {
    DenseMatvec(dfdy, NV_DATA_S(ys[i]), NV_DATA_S(ySdot[i]));
  }

  // the derivative of the rhs wrt the parameter
  assert(data->wrtIn>=0);
  const boost::any& dfdp_any = data->rhs->Jacobian(data->wrtIn, 0, ref_vector<boost::any>(data->inputs.begin(), data->inputs.begin()+data->rhs->numInputs));
  const DlsMat& dfdp = boost::any_cast<const DlsMat&>(dfdp_any);

  // loop through and fill in the rhs vectors stored in ySdot
  N_Vector col = N_VNew_Serial(dfdp->M);
  for( unsigned int i=0; i<Ns; ++i ) {
    NV_DATA_S(col) = DENSE_COL(dfdp, i);

    N_VLinearSum(1.0, ySdot[i], 1.0, col, ySdot[i]);
  }
  NV_DATA_S(col) = nullptr;
  N_VDestroy(col);

  return 0;
}

std::vector<std::pair<unsigned int, unsigned int> > ODEBase::TimeIndices(ref_vector<boost::any> const& outputTimes) {
  // the number of outputs
  const unsigned int nOuts = outputTimes.size();
  
  // each element corresponds to a vector of desired times, first: the current index of that vector, second: the size of that vector
  std::vector<std::pair<unsigned int, unsigned int> > timeIndices(nOuts);

  // resize the output vector
  const boost::any none = boost::none;
  outputs.resize(nOuts, none);

  // loop through the desired outputs
  for( unsigned int i=0; i<nOuts; ++i ) { 
    timeIndices[i].first = 0; // the first index is zero ...
    timeIndices[i].second = algebra->Size(outputTimes[i]); // the size of the vector

    // the the size is of this output >1, set the size of that output vector
    if( timeIndices[i].second>1 ) {
      outputs[i] = std::vector<N_Vector>(timeIndices[i].second);
    } else {
      outputs[i] = N_Vector();
    }
  }

  // return as reference
  return timeIndices;
}

bool ODEBase::NextTime(std::pair<double, int>& nextTime, std::vector<std::pair<unsigned int, unsigned int> >& timeIndices, ref_vector<boost::any> const& outputTimes) const {
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
    const double t = boost::any_cast<double>(algebra->AccessElement(outputTimes[i], timeIndices[i].first));

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

void ODEBase::ClearOutputs() {
  // destory the outputs
  for( unsigned int i=0; i<outputs.size(); ++i ) {
    // check if it is a N_Vector
    if( N_VectorName.compare(outputs[i].type().name())==0 ) {
      N_Vector& vec = boost::any_cast<N_Vector&>(outputs[i]);
      N_VDestroy(vec);
    }
    // check if it is a std::vector<N_Vector>
    if( stdvecN_VectorName.compare(outputs[i].type().name())==0 ) {
      std::vector<N_Vector>& vec = boost::any_cast<std::vector<N_Vector>&>(outputs[i]);
      for( auto it : vec ) {
	N_VDestroy(it);
      }
    }
  }
  
  // clear the outputs
  outputs.clear();
}

void ODEBase::ClearJacobian() {
  // check the jacobian
  if( jacobian ) {
    // check if it is a N_Vector
    if( DlsMatName.compare((*jacobian).type().name())==0 ) {
      DestroyMat(boost::any_cast<DlsMat&>(*jacobian));
    }
    // check if it is a std::vector<DlsMat>
    if( stdvecDlsMatName.compare((*jacobian).type().name())==0 ) {
      std::vector<DlsMat>& vec = boost::any_cast<std::vector<DlsMat>&>((*jacobian));
      for( auto it : vec ) {
	DestroyMat(it);
      }
    }
    
    // reset the jacobian
    jacobian = boost::none;
  }
}

void ODEBase::ClearJacobianAction() {
  // destory the jacobianAction
  if( jacobianAction ) {
    // check if it is a N_Vector
    if( N_VectorName.compare((*jacobianAction).type().name())==0 ) {
      N_VDestroy(boost::any_cast<N_Vector&>(*jacobianAction));
    }
    // check if it is a std::vector<N_Vector>
    if( stdvecN_VectorName.compare((*jacobianAction).type().name())==0 ) {
      std::vector<N_Vector>& vec = boost::any_cast<std::vector<N_Vector>&>((*jacobianAction));
      for( auto it : vec ) {
	N_VDestroy(it);
      }
    }
  }

  // reset the jacobianAction
  jacobianAction = boost::none;
}

void ODEBase::ClearJacobianTransposeAction() {
  // destory the jacobianTransposeAction
  if( jacobianTransposeAction ) {
    // check if it is a N_Vector
    if( N_VectorName.compare((*jacobianTransposeAction).type().name())==0 ) {
      N_VDestroy(boost::any_cast<N_Vector&>(*jacobianTransposeAction));
    }
    // check if it is a std::vector<N_Vector>
    if( stdvecN_VectorName.compare((*jacobianTransposeAction).type().name())==0 ) {
      std::vector<N_Vector>& vec = boost::any_cast<std::vector<N_Vector>&>((*jacobianTransposeAction));
      for( auto it : vec ) {
	N_VDestroy(it);
      }
    }
  }

  // reset the jacobianTransposeAction
  jacobianTransposeAction = boost::none;
}

void ODEBase::ClearResults() {
  ClearOutputs();
  ClearJacobian();
  ClearJacobianAction();
  ClearJacobianTransposeAction();
}
