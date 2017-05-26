#include "MUQ/Modeling/ODE.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;

ODE::ODE(std::shared_ptr<WorkPiece> rhs, pt::ptree const& pt, std::shared_ptr<AnyAlgebra> algebra) : ODEBase(rhs, pt, algebra) {}

ODE::~ODE() {}

void ODE::Integrate(ref_vector<boost::any> const& inputs) const {
  // the number of inputs must be greater than the number of inputs required by the rhs
  assert(inputs.size()>rhs->numInputs);
  
  // get the dimension of the state
  const unsigned int dim = algebra->VectorDimensionBase(inputs[0]);
  assert(dim>0);
  
  // initialize the state (set the initial conditions)
  N_Vector state = nullptr;
  InitializeState(state, inputs[0], dim);

  // crete the solver memory
  void* cvode_mem = nullptr;
  CreateSolverMemory(cvode_mem);

  // each element corresponds to a vector of desired times, first: the current index of that vector, second: the size of that vector
  std::vector<std::pair<unsigned int, unsigned int> > timeIndices(inputs.size()-rhs->numInputs);
  for( unsigned int i=0; i<timeIndices.size(); ++i ) { // loop through the desired outputs
    timeIndices[i].first = 0; // start at the first index
    timeIndices[i].second = algebra->VectorDimensionBase(inputs[rhs->numInputs+i]); // the size of the vector
  }

  // first: the next time to integrate to, second: the output index
  std::pair<double, int> nextTime;

  // the current time
  double tcurrent = 0.0;

  // while we should keep integrating
  while( NextTime(nextTime, timeIndices, ref_vector<boost::any>(inputs.begin()+rhs->numInputs, inputs.end())) ) {
    if( fabs(nextTime.first-tcurrent)>1.0e-14 ) { // we have to move forward in time
      std::cout << nextTime.first << " " << nextTime.second << std::endl;
    } else {
    }
  }
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
