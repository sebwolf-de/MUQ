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

ODEBase::ODEBase(std::shared_ptr<WorkPiece> rhs, pt::ptree const& pt, std::shared_ptr<AnyAlgebra> algebra) : WorkPiece(), rhs(rhs), algebra(algebra), linSolver(pt.get<std::string>("ODESolver.LinearSolver", "Dense")), reltol(pt.get<double>("ODESolver.RelativeTolerance", 1.0e-8)), abstol(pt.get<double>("ODESolver.AbsoluteTolerance", 1.0e-8)), maxStepSize(pt.get<double>("ODESolver.MaxStepSize", 1.0)) {
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
  
  // the type of the first input (the state) for the rhs
  assert(stateType.compare(rhs->InputType(0, false))==0);

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
    NV_Ith_S(state, i) = boost::any_cast<double>(algebra->AccessElementBase(i, ic));
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
  data->inputs[0] = anyref;

  // evaluate the rhs
  const std::vector<boost::any>& result = data->rhs->Evaluate(ref_vector<boost::any>(data->inputs.begin(), data->inputs.begin()+data->rhs->numInputs));
  NV_DATA_S(statedot) = NV_DATA_S(boost::any_cast<N_Vector>(result[0]));

  return 0;
}

int ODEBase::RHSJacobianAction(N_Vector v, N_Vector Jv, realtype time, N_Vector state, N_Vector rhs, void *user_data, N_Vector tmp) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);

  // set the state input
  const boost::any& anyref = state;
  data->inputs[0] = anyref;

  // compute the jacobain wrt the state
  const boost::any& jacobianAction = data->rhs->JacobianAction(0, 0, v, ref_vector<boost::any>(data->inputs.begin(), data->inputs.begin()+data->rhs->numInputs));
  NV_DATA_S(Jv) = NV_DATA_S(boost::any_cast<N_Vector>(jacobianAction));

  return 0;
}

int ODEBase::RHSJacobian(long int N, realtype time, N_Vector state, N_Vector rhs, DlsMat jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  // get the data type
  ODEData* data = (ODEData*)user_data;
  assert(data);
  assert(data->rhs);
  
  // set the state input
  const boost::any& anyref = state;
  data->inputs[0] = anyref;

  // evaluate the jacobian
  const boost::any& jcbn = data->rhs->Jacobian(0, 0, ref_vector<boost::any>(data->inputs.begin(), data->inputs.begin()+data->rhs->numInputs));
  jac = boost::any_cast<DlsMat>(jcbn);

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

  // the derivative of the rhs wrt the parameter
  assert(data->wrtIn>=0);
  const boost::any& dfdp_any = data->rhs->Jacobian(data->wrtIn, 0, ref_vector<boost::any>(data->inputs.begin(), data->inputs.begin()+data->rhs->numInputs));
  const DlsMat& dfdp = boost::any_cast<const DlsMat&>(dfdp_any);

  // loop through and fill in the rhs vectors stored in ySdot
  for( unsigned int i=0; i<Ns; ++i ) {
    N_Vector vec = N_VNew_Serial(dfdy->M);
    DenseMatvec(dfdy, NV_DATA_S(ys[i]), NV_DATA_S(vec));

    N_Vector col = N_VNew_Serial(dfdp->M);
    NV_DATA_S(col) = DENSE_COL(dfdp, i);

    N_VLinearSum(1.0, vec, 1.0, col, ySdot[i]);
  }

  return 0;
}
