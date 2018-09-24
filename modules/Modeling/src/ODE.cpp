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

ODE::ODE(std::shared_ptr<ModPiece> rhs, pt::ptree const& pt) : ODEBase(rhs, InputSizes(rhs, pt), OutputSizes(rhs, pt), pt) {}

ODE::~ODE() {}

Eigen::VectorXi ODE::InputSizes(std::shared_ptr<ModPiece> rhs,  boost::property_tree::ptree const& pt) {
  // number of inputs for the RHS plus the time vector
  Eigen::VectorXi inSizes(rhs->numInputs+1);

  // same as rhs
  inSizes.head(rhs->numInputs) = rhs->inputSizes;

  // number of times
  inSizes(rhs->numInputs) = pt.get<unsigned int>("NumObservations");

  return inSizes;
}

Eigen::VectorXi ODE::OutputSizes(std::shared_ptr<ModPiece> rhs,  boost::property_tree::ptree const& pt) {
  // each output time is an output
  const unsigned int numOuts = pt.get<unsigned int>("NumObservations");

  // the index of the state for the rhs is 0 if we are autonomous and 1 if we are not
  const unsigned int rhsStateIndex = pt.get<bool>("Autonomous", true) ? 0 : 1;

  // number of elements for each output is the state size
  const unsigned int stateSize = rhs->inputSizes(rhsStateIndex);

  // number of inputs for the RHS plus the time vector
  return Eigen::VectorXi::Constant(1, numOuts*stateSize);
}


void ODE::Integrate(ref_vector<Eigen::VectorXd> const& inputs, int const wrtIn, N_Vector const& vec, DerivativeMode const& mode) {
  // the number of inputs must be greater than the number of inputs required by the rhs
  assert(inputs.size()>rhs->numInputs-(autonomous? 0 : 1));

  // get the state size
  const unsigned int stateSize = autonomous? inputSizes(0) : inputSizes(1);

  // create the state vector
  N_Vector state = N_VNew_Serial(stateSize);
  Eigen::Map<Eigen::VectorXd> stateMap(NV_DATA_S(state), stateSize);
  stateMap = inputs.at(0).get();

  // create a data structure to pass around in Sundials
  auto data = std::make_shared<ODEData>(rhs, ref_vector<Eigen::VectorXd>(inputs.begin(), inputs.begin()+rhs->numInputs), autonomous, wrtIn);
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

  if( wrtIn>=0 ) { // we are computing the forward sensitivity
    // integrate forward in time with dervative information
    ForwardSensitivity(cvode_mem, state, wrtIn, inputs.back(), ref_vector<Eigen::VectorXd>(inputs.begin(), inputs.begin()+rhs->numInputs));

    // free integrator memory
    CVodeFree(&cvode_mem);
    N_VDestroy(state);

    return;
  }

  // integrate forward in time
  ForwardIntegration(cvode_mem, state, inputs.back());

  // free integrator memory
  CVodeFree(&cvode_mem);
  N_VDestroy(state);
}

void ODE::ForwardIntegration(void *cvode_mem, N_Vector& state, Eigen::VectorXd const& outputTimes) {
  outputs.resize(1);
  outputs[0] = Eigen::VectorXd::Constant(outputSizes(0), std::numeric_limits<double>::quiet_NaN());

  // start the clock at t=0
  double tcurrent = 0.0;

  // loop through each time
  for( unsigned int t=0; t<outputTimes.size(); ++t ) {
    // if we have to move forward --- i.e., not at the initial time or another output does not need the current time
    if( std::fabs(outputTimes(t)-tcurrent)>1.0e-14 ) { // we have to move forward in time
      int flag = CVode(cvode_mem, outputTimes(t), state, &tcurrent, CV_NORMAL);
      assert(CheckFlag(&flag, "CVode", 1));
    }

    Eigen::Map<Eigen::VectorXd> stateref(NV_DATA_S(state), NV_LENGTH_S(state));
    outputs[0].segment(t*inputSizes(0), inputSizes(0)) = stateref;
  }
}

void ODE::ForwardSensitivity(void *cvode_mem, N_Vector& state, unsigned int const wrtIn, Eigen::VectorXd const& outputTimes, ref_vector<Eigen::VectorXd> const& rhsInputs) {
  // compute the parameter size
  unsigned int paramSize = inputSizes(wrtIn);

  // initialize the jacobian to zero
  jacobian = Eigen::MatrixXd::Zero(outputSizes(0), paramSize);

  // sensState will hold several N_Vectors (because it's a Jacobian)
  N_Vector *sensState = nullptr;

  // set up sensitivity vector
  sensState = N_VCloneVectorArray_Serial(paramSize, state);
  assert(CheckFlag((void *)sensState, "N_VCloneVectorArray_Serial", 0));

  // initialize the sensitivies to zero
  for( int is=0; is<paramSize; ++is ) {
    N_VConst(0.0, sensState[is]);
  }

  // set up solver for sensitivity
  SetUpSensitivity(cvode_mem, paramSize, sensState);

  // start the clock at t=0
  double tcurrent = 0.0;

  // loop through each time
  for( unsigned int t=0; t<outputTimes.size(); ++t ) {
    jacobian.block(t*inputSizes(0), 0, inputSizes(0), paramSize) = wrtIn==0? (Eigen::MatrixXd)Eigen::MatrixXd::Identity(NV_LENGTH_S(state), paramSize) : (Eigen::MatrixXd)Eigen::MatrixXd::Zero(NV_LENGTH_S(state), paramSize);

    // if we have to move forward --- i.e., not at the initial time
    if( std::fabs(outputTimes(t)-tcurrent)>1.0e-14 ) {
      // integrate until we hit the next time
      int flag = CVode(cvode_mem, outputTimes(t), state, &tcurrent, CV_NORMAL);
      assert(CheckFlag(&flag, "CVode", 1));

      // get the sensitivities if we are differentating wrt an input to the rhs
      if( wrtIn!=numInputs-1 ) {
        // get the sensitivities
        int flag = CVodeGetSens(cvode_mem, &tcurrent, sensState);
        assert(CheckFlag(&flag, "CVodeGetSense", 1));

        for( unsigned int i=0; i<NV_LENGTH_S(state); ++i ) {
          for( unsigned int j=0; j<paramSize; ++j ) {
            jacobian(t*inputSizes(0)+i, j) += NV_Ith_S(sensState[j], i);
          }
        }
      }
    }
  }

  // destroy the sensivity
  if( sensState ) {
    N_VDestroyVectorArray_Serial(sensState, paramSize);
  }
}

void ODE::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  // integrate forward in time
  Integrate(inputs);
}

void ODE::JacobianImpl(unsigned int const wrtOut, unsigned int const wrtIn, ref_vector<Eigen::VectorXd> const& inputs) {
  assert(wrtOut==0);
  if( wrtIn!=numInputs-1 ) {
    // integrate forward in time, using forward senstivities
    return Integrate(inputs, wrtIn);
  }

  Evaluate(inputs);

  // initialize the jacobian to zero
  jacobian = Eigen::MatrixXd::Zero(outputSizes(0), inputSizes(wrtIn));

  ref_vector<Eigen::VectorXd> ins(inputs.begin(), inputs.end()-1);
  Eigen::VectorXd time = Eigen::VectorXd::Zero(1);
  if( !autonomous ) {
    ins.insert(ins.begin(), std::cref(time));
  }

  // loop through each time
  for( unsigned int t=0; t<inputs.back().get().size(); ++t ) {
    const Eigen::VectorXd& state = outputs[0].segment(inputSizes(0)*t, inputSizes(0));

    if( autonomous ) {
      ins[0] = std::cref(state);
    } else {
      time(0) = inputs.back() (t);
      ins[0] = std::cref(time);
      ins[1] = std::cref(state);
    }

    // evaluate the right hand side
    const Eigen::VectorXd& eval = rhs->Evaluate(ins) [0];
    for( unsigned int i=0; i<inputSizes(0); ++i ) {
      jacobian(inputSizes(0)*t+i, t) = eval(i);
    }
  }
}

void ODE::GradientImpl(unsigned int const wrtOut, unsigned int const wrtIn, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& vec) {
  assert(wrtOut==0);
  if( wrtIn!=numInputs-1 ) { return BackwardIntegration(inputs, wrtIn, vec); }

  Evaluate(inputs);

  // initialize the gradient to zero
  gradient = Eigen::VectorXd::Zero(inputSizes(wrtIn));

  ref_vector<Eigen::VectorXd> ins(inputs.begin(), inputs.end()-1);
  Eigen::VectorXd time = Eigen::VectorXd::Zero(1);
  if( !autonomous ) {
    ins.insert(ins.begin(), std::cref(time));
  }

  // loop through each time
  for( unsigned int t=0; t<inputs.back().get().size(); ++t ) {
    const Eigen::VectorXd& state = outputs[0].segment(inputSizes(0)*t, inputSizes(0));

    if( autonomous ) {
      ins[0] = std::cref(state);
    } else {
      time(0) = inputs.back() (t);
      ins[0] = std::cref(time);
      ins[1] = std::cref(state);
    }

    // evaluate the right hand side
    const Eigen::VectorXd& eval = rhs->Evaluate(ins) [0];
    gradient(t) = eval.dot(vec.segment(t*inputSizes(0), inputSizes(0)));
  }
}

void ODE::BackwardIntegration(ref_vector<Eigen::VectorXd> const& inputs, unsigned int const wrtIn, Eigen::VectorXd const& vec) {
  // get the output times
  const Eigen::VectorXd& outputTimes = inputs.back();
  double finalTime = outputTimes[outputTimes.size()-1];

  // get the sizes
  const unsigned int stateSize = inputSizes(0);
  const unsigned int paramSize = inputSizes(wrtIn);

  // create the state vector
  N_Vector state = N_VNew_Serial(stateSize);
  Eigen::Map<Eigen::VectorXd> stateMap(NV_DATA_S(state), stateSize);
  stateMap = inputs.at(0).get();

  // create the adjoint vector
  N_Vector lambda = N_VNew_Serial(stateSize);
  N_Vector nvGrad = N_VNew_Serial(paramSize);
  Eigen::Map<Eigen::VectorXd> lambdaMap(NV_DATA_S(lambda), stateSize);
  Eigen::Map<Eigen::VectorXd> nvGradMap(NV_DATA_S(nvGrad), paramSize);
  lambdaMap = -1.0*vec.tail(stateSize);
  nvGradMap = Eigen::VectorXd::Zero(paramSize);

  // create a data structure to pass around in Sundials
  auto data = std::make_shared<ODEData>(rhs, ref_vector<Eigen::VectorXd>(inputs.begin(), inputs.begin()+rhs->numInputs), autonomous, wrtIn);
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
  const int indexB = CreateSolverMemoryB(cvode_mem, finalTime, lambda, nvGrad, data);

  // start the clock at t=0
  double tcurrent = 0.0;
  int numCheckpts;

  // forward integration
  for( unsigned int t=0; t<outputTimes.size(); ++t ) {
    if( std::fabs(outputTimes(t)-tcurrent)>1.0e-14 ) {
      int flag = CVodeF(cvode_mem, outputTimes(t), state, &tcurrent, CV_NORMAL, &numCheckpts);
      assert(CheckFlag(&flag, "CVodeF", 1));
    }
  }

  for( int t=outputTimes.size()-2; t>=0; --t ) {
    if( std::fabs(outputTimes(t)-tcurrent)>1.0e-14 ) {
      // integrate the adjoint variable back in time
      int flag = CVodeB(cvode_mem, outputTimes(t), CV_NORMAL);
      assert(CheckFlag(&flag, "CVodeB", 1));

      // get the adjoint variable
      flag = CVodeGetB(cvode_mem, indexB, &tcurrent, lambda);
      assert(CheckFlag(&flag, "CVodeGetB", 1));
    }

    // get the quadrature
    int flag = CVodeGetQuadB(cvode_mem, indexB, &tcurrent, nvGrad);
    assert(CheckFlag(&flag, "CVodeGetQuadB", 1));

    // discontinuous change in variable
    lambdaMap -= vec.segment(t*stateSize, stateSize);

    // reinitialize the adjoint integrator because of the discontinuity in the state
    flag = CVodeReInitB(cvode_mem, indexB, outputTimes(t), lambda);
    assert(CheckFlag(&flag, "CVodeReInitB", 1));
    flag = CVodeQuadReInitB(cvode_mem, indexB, nvGrad);
    assert(CheckFlag(&flag, "CVodeQuadReInitB", 1));
  }

  // makes sure we integrate all the way back to zero
  if( std::fabs(tcurrent)>1.0e-14 ) {
    // integrate the adjoint variable back in time
    int flag = CVodeB(cvode_mem, 0.0, CV_NORMAL);
    assert(CheckFlag(&flag, "CVodeB", 1));

    // get the adjoint variable
    flag = CVodeGetB(cvode_mem, indexB, &tcurrent, lambda);
    assert(CheckFlag(&flag, "CVodeGetB", 1));
  }

  // get the gradient
  int flag = CVodeGetQuadB(cvode_mem, indexB, &tcurrent, nvGrad);
  assert(CheckFlag(&flag, "CVodeGetQuadB", 1));
  gradient = nvGradMap;

  // if wrt to the initial conditions add identity*sens
  if( wrtIn==0 ) {
    for( unsigned int t=0; t<outputTimes.size(); ++t ) {
      gradient += vec.segment(t*stateSize, stateSize);
    }
  }

  // free memory
  N_VDestroy_Serial(state);
  N_VDestroy_Serial(lambda);
  N_VDestroy_Serial(nvGrad);

  CVodeFree(&cvode_mem);
}
