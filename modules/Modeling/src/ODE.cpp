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

ODE::ODE(std::shared_ptr<ModPiece> rhs, pt::ptree const& pt) : ODEBase(rhs, pt) {}

ODE::~ODE() {}

void ODE::Integrate(ref_vector<Eigen::VectorXd> const& inputs, int const wrtIn, int const wrtOut, N_Vector const& vec, DerivativeMode const& mode) {
  // the number of inputs must be greater than the number of inputs required by the rhs
  assert(inputs.size()>rhs->numInputs-(autonomous? 0 : 1));

  // clear the results
  ClearResults();

  // create the state vector (have to do a hard copy --- N_Vector is a pointer to the data, the pointer has been declared const, not the data)
  N_Vector state;
  DeepCopy(state, boost::any_cast<const N_Vector&>(inputs[0]));

  // create a data structure to pass around in Sundials
  auto data = std::make_shared<ODEData>(rhs, inputs, autonomous, wrtIn, wrtOut);
  if( !autonomous ) {
    data->inputs.insert(data->inputs.begin(), Eigen::VectorXd::Zero(1));
  }

  // set solver to null
  void* cvode_mem = nullptr;

  // create the solver memory
  cvode_mem = CVodeCreate(multiStep, solveMethod);
  assert(CheckFlag((void*)cvode_mem, "CVodeCreate", 0));

  // initialize the solver
  CreateSolverMemory(cvode_mem, state, data);

  if( wrtIn>=0 && wrtOut>=0 ) { // we are computing the forward sensitivity
    // compute the parameter size
    unsigned int paramSize = 1;
    if( wrtIn<rhs->numInputs ) {
      paramSize = inputs[wrtIn].get().size();
    }

    // integrate forward in time with dervative information
    ForwardSensitivity(cvode_mem, state, paramSize, wrtIn, inputs[rhs->numInputs+wrtOut], ref_vector<Eigen::VectorXd>(inputs.begin(), inputs.begin()+rhs->numInputs), mode, vec);

    // free integrator memory
    CVodeFree(&cvode_mem);
    N_VDestroy(state);

    return;
  }

  // integrate forward in time
  ForwardIntegration(cvode_mem, state, ref_vector<Eigen::VectorXd>(inputs.begin()+rhs->numInputs, inputs.end()));

  // free integrator memory
  CVodeFree(&cvode_mem);
  N_VDestroy(state);
}

void ODE::ForwardIntegration(void *cvode_mem, N_Vector& state, ref_vector<Eigen::VectorXd> const& outputTimes) {
  /*// each element corresponds to a vector of desired times, first: the current index of that vector, second: the size of that vector
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
  }*/
}

void ODE::ForwardSensitivity(void *cvode_mem, N_Vector& state, unsigned int const paramSize, unsigned int const wrtIn, boost::any const& outputTimes, ref_vector<Eigen::VectorXd> const& rhsInputs, DerivativeMode const& mode, N_Vector const& vec) {
  /*// sensState will hold several N_Vectors (because it's a Jacobian)
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
  const unsigned int ntimes = algebra->Size(outputTimes);

  // initialize the derivative information (jacobian jacobianAction, jacobianTransposeAction)
  InitializeDerivative(ntimes, NV_LENGTH_S(state), paramSize, mode);

  // start the clock at t=0
  double tcurrent = 0.0;

  // loop through each time
  for( unsigned int i=0; i<ntimes; ++i ) {
    // the next time
    const double time = boost::any_cast<double>(algebra->AccessElement(outputTimes, i));

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

  // destroy the sensivity
  if( sensState ) {
    N_VDestroyVectorArray_Serial(sensState, paramSize);
  }*/
}

void ODE::SaveDerivative(unsigned int const ntimes, unsigned int const timeIndex, unsigned int const paramSize, unsigned int const wrtIn, N_Vector* sensState, N_Vector const& state, ref_vector<boost::any> rhsInputs, N_Vector const& vec, DerivativeMode const& mode) {
  /*if( ntimes==1 ) {
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
    case DerivativeMode::Jac: { // save the jacobian
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
  }*/
}

void ODE::SaveJacobianAction(N_Vector& jacAct, unsigned int const ncols, unsigned int const wrtIn, N_Vector* sensState, N_Vector const& state, ref_vector<boost::any> rhsInputs, N_Vector const& vec) const {
  /*// the number of rows in the Jacobian
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
  }*/
}

void ODE::SaveJacobianTransposeAction(N_Vector& jacTransAct, unsigned int const ncols, unsigned int const wrtIn, N_Vector* sensState, N_Vector const& state, ref_vector<boost::any> rhsInputs, N_Vector const& vec) const {
  /*// the number of rows in the Jacobian
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
  }*/
}

void ODE::SaveJacobian(DlsMat& jac, unsigned int const ncols, unsigned int const wrtIn, N_Vector *sensState, N_Vector const& state, ref_vector<boost::any> rhsInputs) const {
  /*// the number of rows in the Jacobian
  const unsigned int nrows = NV_LENGTH_S(state);

  if( wrtIn>=rhs->numInputs ) {
    // set the state input
    const boost::any& anyref = state;
    rhsInputs[0] = anyref;

    // evaluate the right hand side
    const std::vector<boost::any>& eval = rhs->Evaluate(rhsInputs);

    // save as a matrix (need a deep copy since eval will be destory next time rhs->Evaluate is called)
    jac = NewDenseMat(nrows, 1);
    const N_Vector& val = boost::any_cast<const N_Vector&>(eval[0]);
    for( unsigned int i=0; i<nrows; ++i ) {
      DENSE_ELEM(jac, i, 0) = NV_Ith_S(val, i);
    }

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
  }*/
}

void ODE::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  // integrate forward in time
  Integrate(inputs);
}

void ODE::JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<Eigen::VectorXd> const& inputs) {
  // integrate forward in time, using forward senstivities
  Integrate(inputs, wrtIn, wrtOut);
}

void ODE::ApplyJacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& vec) {
    /*// check the size of the vector
  if( wrtIn<rhs->numInputs ) {
    assert(NV_LENGTH_S(boost::any_cast<const N_Vector&>(vec))==algebra->Size(inputs[wrtIn]));
  } else {
    assert(NV_LENGTH_S(boost::any_cast<const N_Vector&>(vec))==1);
  }

  // integrate forward in time, using forward senstivities
  Integrate(inputs, wrtIn, wrtOut, boost::any_cast<N_Vector>(vec), DerivativeMode::JacAction);*/
}

void ODE::GradientImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& vec) {
  // check the size of the vector
  assert(NV_LENGTH_S(boost::any_cast<const N_Vector&>(vec))==NV_LENGTH_S(boost::any_cast<const N_Vector&>(inputs[0])));

  // integrate forward in time, using forward senstivities
  Integrate(inputs, wrtIn, wrtOut, boost::any_cast<N_Vector>(vec), DerivativeMode::JacTransAction);
}
