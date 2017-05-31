#include "MUQ/Modeling/ODE.h"

#include <cvodes/cvodes.h> // prototypes for CVODE fcts. and consts. 
#include <cvodes/cvodes_spgmr.h> // prototypes & constants for CVSPGMR solver 
#include <cvodes/cvodes_spbcgs.h> // prototypes & constants for CVSPBCG solver 
#include <cvodes/cvodes_sptfqmr.h> // prototypes & constants for SPTFQMR solver 
#include <cvodes/cvodes_dense.h> // prototype for CVDense 

#include <sundials/sundials_dense.h> /* use generic DENSE solver in preconditioning */
#include <sundials/sundials_types.h> // definition of type 
#include <sundials/sundials_math.h>  // contains the macros ABS, SQR, and EXP 
#include <sundials/sundials_direct.h>

namespace pt = boost::property_tree;
using namespace muq::Modeling;

ODE::ODE(std::shared_ptr<WorkPiece> rhs, pt::ptree const& pt, std::shared_ptr<AnyAlgebra> algebra) : ODEBase(rhs, pt, algebra) {}

ODE::~ODE() {}

void ODE::Integrate(ref_vector<boost::any> const& inputs, int const wrtIn, int const wrtOut, N_Vector const& vec, DerivativeMode const& mode) {
  // the number of inputs must be greater than the number of inputs required by the rhs
  assert(inputs.size()>rhs->numInputs);

  // create the state vector (have to do a hard copy --- N_Vector is a pointer to the data, the pointer has been declared const, not the data)
  N_Vector state;
  DeepCopy(state, boost::any_cast<const N_Vector&>(inputs[0]));

  // create a data structure to pass around in Sundials
  auto data = std::make_shared<ODEData>(rhs, inputs, wrtIn, wrtOut);

  // set solver to null
  void* cvode_mem = nullptr;

  // create the solver memory
  cvode_mem = CVodeCreate(multiStep, solveMethod);
  assert(CheckFlag((void*)cvode_mem, "CVodeCreate", 0));

  // initialize the solver
  CreateSolverMemory(cvode_mem, state, data);
  
  if( wrtIn>=0 && wrtOut>=0 ) { // we are computing the forward sensitivity
    // comptue the parameter size
    unsigned int paramSize = 1;
    if( wrtIn<rhs->numInputs ) {
      paramSize = algebra->VectorDimensionBase(inputs[wrtIn]);
    }

    // integrate forward in time with dervative information
    ForwardSensitivity(cvode_mem, state, paramSize, wrtIn, inputs[rhs->numInputs+wrtOut], ref_vector<boost::any>(inputs.begin(), inputs.begin()+rhs->numInputs), mode, vec);

    return;
  }

  // integrate forward in time
  ForwardIntegration(cvode_mem, state, ref_vector<boost::any>(inputs.begin()+rhs->numInputs, inputs.end()));
}

void ODE::ForwardIntegration(void *cvode_mem, N_Vector& state, ref_vector<boost::any> const& outputTimes) {
  // each element corresponds to a vector of desired times, first: the current index of that vector, second: the size of that vector
  std::vector<std::pair<unsigned int, unsigned int> > timeIndices = TimeIndices(outputTimes);

  // first: the next time to integrate to, second: the output index
  std::pair<double, int> nextTime;

  // start the clock at t=0
  double tcurrent = 0.0;

  // loop through each time 
  while( NextTime(nextTime, timeIndices, outputTimes) ) {
    // if we have to move forward --- i.e., not at the initial time or another output does not need the current time
    if( std::fabs(nextTime.first-tcurrent)>1.0e-14 ) { // we have to move forward in time
      int flag = CVode(cvode_mem, nextTime.first, state, &tcurrent, CV_NORMAL);
      assert(CheckFlag(&flag, "CVode", 1));
    }

    // save the result at this timestep
    if( timeIndices[nextTime.second].second>1 ) { // the output has more than one compnent ...
      // .. save the current state as an element in the vector
      DeepCopy(boost::any_cast<std::vector<N_Vector>&>(outputs[nextTime.second]) [timeIndices[nextTime.second].first-1], state);
    } else { // the out has one component ...
      // ... save the current state (not inside a vector)
      DeepCopy(boost::any_cast<N_Vector&>(outputs[nextTime.second]), state);
    }
  }
}

void ODE::ForwardSensitivity(void *cvode_mem, N_Vector& state, unsigned int const paramSize, unsigned int const wrtIn, boost::any const& outputTimes, ref_vector<boost::any> const& rhsInputs, DerivativeMode const& mode, N_Vector const& vec) {
  // sensState will hold several N_Vectors (because it's a Jacobian)
  N_Vector *sensState = nullptr;

  // if we are computing the derivative wrt an input to the rhs ...
  if( wrtIn<rhs->numInputs ) {
    // set up sensitivity vector
    sensState = N_VCloneVectorArray_Serial(paramSize, state);
    assert(CheckFlag((void *)sensState, "N_VCloneVectorArray_Serial", 0));
    
    // initialize the sensitivies to zero
    for( int is=0; is<paramSize; ++is ) {
      N_VConst(0.0, sensState[is]);
    }
    
    // set up solver for sensitivity
    SetUpSensitivity(cvode_mem, paramSize, sensState);
  }
  
  // number of output times
  const unsigned int ntimes = algebra->VectorDimensionBase(outputTimes);

  // initialize the derivative information (jacobian jacobianAction, jacobianTransposeAction)
  InitializeDerivative(ntimes, NV_LENGTH_S(state), paramSize, mode);

  // start the clock at t=0
  double tcurrent = 0.0;

  // loop through each time
  for( unsigned int i=0; i<ntimes; ++i ) {
    // the next time 
    const double time = boost::any_cast<double>(algebra->AccessElementBase(i, outputTimes));

    // if we have to move forward --- i.e., not at the initial time
    if( std::fabs(time-tcurrent)>1.0e-14 ) {
      // integrate until we hit the next time
      int flag = CVode(cvode_mem, time, state, &tcurrent, CV_NORMAL);
      assert(CheckFlag(&flag, "CVode", 1));

      // get the sensitivities if we are differentating wrt an input to the rhs
      if( wrtIn<rhs->numInputs ) {
	// get the sensitivities
	flag = CVodeGetSens(cvode_mem, &tcurrent, sensState);
	assert(CheckFlag(&flag, "CVodeGetSense", 1));
      }
    }

    // save the derivative information
    SaveDerivative(ntimes, i, paramSize, wrtIn, sensState, state, rhsInputs, vec, mode);
  }
}

void ODE::SaveDerivative(unsigned int const ntimes, unsigned int const timeIndex, unsigned int const paramSize, unsigned int const wrtIn, N_Vector* sensState, N_Vector const& state, ref_vector<boost::any> rhsInputs, N_Vector const& vec, DerivativeMode const& mode) {
  if( ntimes==1 ) {
    switch( mode ) { // one output --- matrix or vector
    case DerivativeMode::Jac: { // save the jacobian
      SaveJacobian(boost::any_cast<DlsMat&>(*jacobian), paramSize, wrtIn, sensState, state, rhsInputs);
      return;
    }
    case DerivativeMode::JacAction: { // save the action of the jacobian
      SaveJacobianAction(boost::any_cast<N_Vector&>(*jacobianAction), paramSize, wrtIn, sensState, state, rhsInputs, vec);
      return;
    }
    case DerivativeMode::JacTransAction: { // save the action of the jacobian transpose
      SaveJacobianTransposeAction(boost::any_cast<N_Vector&>(*jacobianTransposeAction), paramSize, wrtIn, sensState, state, rhsInputs, vec);
      return;
    }
    default:
      assert(false);
    }
  } else { // multiple outputs --- vector of matrices or vectors
    switch( mode ) {
    case DerivativeMode::Jac: { // savethe jacobian
     SaveJacobian(boost::any_cast<std::vector<DlsMat>&>(*jacobian) [timeIndex], paramSize, wrtIn, sensState, state, rhsInputs);
     return;
    }
    case DerivativeMode::JacAction: { // save the action of the jacobian
      SaveJacobianAction(boost::any_cast<std::vector<N_Vector>&>(*jacobianAction) [timeIndex], paramSize, wrtIn, sensState, state, rhsInputs, vec);
      return;
    }
    case DerivativeMode::JacTransAction: { // save the action of the jacobian transpose
      SaveJacobianTransposeAction(boost::any_cast<std::vector<N_Vector>&>(*jacobianTransposeAction) [timeIndex], paramSize, wrtIn, sensState, state, rhsInputs, vec);
      return;
    }
    default:
      assert(false);
    }
  }
}

void ODE::SaveJacobianAction(N_Vector& jacAct, unsigned int const ncols, unsigned int const wrtIn, N_Vector* sensState, N_Vector const& state, ref_vector<boost::any> rhsInputs, N_Vector const& vec) const {
  // the number of rows in the Jacobian
  const unsigned int nrows = NV_LENGTH_S(state);

  // create a new jacobian action vector
  jacAct = N_VNew_Serial(nrows);

  if( wrtIn>=rhs->numInputs ) {
    // set the state input
    const boost::any& anyref = state;
    rhsInputs[0] = anyref;

    // evaluate the right had side
    const std::vector<boost::any>& eval = rhs->Evaluate(rhsInputs);

    // multiply by the input vector
    N_VScale(NV_Ith_S(vec, 0), boost::any_cast<const N_Vector&>(eval[0]), jacAct);

    return; 
  }

  // populate the jacobian action
  for( unsigned int row=0; row<nrows; ++row ) {
    NV_Ith_S(jacAct, row) = 0.0;
    for( unsigned int col=0; col<ncols; ++col ) {
      NV_Ith_S(jacAct, row) += NV_Ith_S(sensState[col], row)*NV_Ith_S(vec, col);
    }
  }

  // if it wrt the initial conditions, add the identity
  if( wrtIn==0 ) {
    N_VLinearSum(1.0, jacAct, 1.0, vec, jacAct);
  }
}

void ODE::SaveJacobianTransposeAction(N_Vector& jacTransAct, unsigned int const ncols, unsigned int const wrtIn, N_Vector* sensState, N_Vector const& state, ref_vector<boost::any> rhsInputs, N_Vector const& vec) const {
  // the number of rows in the Jacobian
  const unsigned int nrows = NV_LENGTH_S(state);

  // create a new jacobian action vector
  jacTransAct = N_VNew_Serial(ncols);

  if( wrtIn>=rhs->numInputs ) {
    // set the state input
    const boost::any& anyref = state;
    rhsInputs[0] = anyref;

    // evaluate the right had side
    const std::vector<boost::any>& eval = rhs->Evaluate(rhsInputs);

    // multiply by the input vector
    NV_Ith_S(jacTransAct, 0) = N_VDotProd(vec, boost::any_cast<const N_Vector&>(eval[0]));

    return; 
  }

  // populate the jacobian action
  for( unsigned int col=0; col<ncols; ++col ) {
    NV_Ith_S(jacTransAct, col) = N_VDotProd(vec, sensState[col]);
  }

  // if it wrt the initial conditions, add the identity
  if( wrtIn==0 ) {
    N_VLinearSum(1.0, jacTransAct, 1.0, vec, jacTransAct);
  }
}

void ODE::SaveJacobian(DlsMat& jac, unsigned int const ncols, unsigned int const wrtIn, N_Vector *sensState, N_Vector const& state, ref_vector<boost::any> rhsInputs) const {
  // the number of rows in the Jacobian
  const unsigned int nrows = NV_LENGTH_S(state);
  
  if( wrtIn>=rhs->numInputs ) {
    // set the state input
    const boost::any& anyref = state;
    rhsInputs[0] = anyref;

    // evaluate the right hand side
    const std::vector<boost::any>& eval = rhs->Evaluate(rhsInputs);

    // save as a matrix
    jac = NewDenseMat(nrows, 1);
    DENSE_COL(jac, 0) = NV_DATA_S(boost::any_cast<const N_Vector&>(eval[0]));

    return; 
  }
  
  // create a new jacobian
  jac = NewDenseMat(nrows, ncols);

  // populate the jacobian
  for( unsigned int col=0; col<ncols; ++col ) {
    // need to do a deep copy
    for( unsigned int row=0; row<nrows; ++row ) {
      DENSE_ELEM(jac, row, col) = NV_Ith_S(sensState[col], row);
    }
  }

  // if it wrt the initial conditions, add the identity
  if( wrtIn==0 ) {
    // the matrix is square
    assert(jac->M==jac->N);
    AddIdentity(jac);
  }
}

void ODE::SetUpSensitivity(void *cvode_mem, unsigned int const paramSize, N_Vector *sensState) const {
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

void ODE::InitializeDerivative(unsigned int const ntimes, unsigned int const stateSize, unsigned int const paramSize, DerivativeMode const& mode) {
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

std::vector<std::pair<unsigned int, unsigned int> > ODE::TimeIndices(ref_vector<boost::any> const& outputTimes) {
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
    timeIndices[i].second = algebra->VectorDimensionBase(outputTimes[i]); // the size of the vector

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
  // integrate forward in time
  Integrate(inputs);
}

void ODE::JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) {
  // integrate forward in time, using forward senstivities
  Integrate(inputs, wrtIn, wrtOut);
}

void ODE::JacobianActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) {
  // check the size of the vector
  if( wrtIn<rhs->numInputs ) {
    assert(NV_LENGTH_S(boost::any_cast<const N_Vector&>(vec))==algebra->VectorDimensionBase(inputs[wrtIn]));
  } else {
    assert(NV_LENGTH_S(boost::any_cast<const N_Vector&>(vec))==1);
  }

  // integrate forward in time, using forward senstivities
  Integrate(inputs, wrtIn, wrtOut, boost::any_cast<N_Vector>(vec), DerivativeMode::JacAction);
}

void ODE::JacobianTransposeActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) {
  // check the size of the vector
  assert(NV_LENGTH_S(boost::any_cast<const N_Vector&>(vec))==NV_LENGTH_S(boost::any_cast<const N_Vector&>(inputs[0])));

  // integrate forward in time, using forward senstivities
  Integrate(inputs, wrtIn, wrtOut, boost::any_cast<N_Vector>(vec), DerivativeMode::JacTransAction);
}

int ODE::ForwardSensitivityRHS(int Ns, realtype time, N_Vector y, N_Vector ydot, N_Vector *ys, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2) {
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
