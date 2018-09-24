#include "MUQ/Modeling/ODEData.h"

#include "boost/none.hpp"

using namespace muq::Modeling;

// construct basic ode data
ODEData::ODEData(std::shared_ptr<ModPiece> rhs, ref_vector<Eigen::VectorXd> const& refinputs, bool const autonomous, int const wrtIn) : rhs(rhs), autonomous(autonomous), wrtIn(wrtIn) {
  inputs.resize(refinputs.size());
  for( unsigned int i=0; i<refinputs.size(); ++i ) { inputs[i] = refinputs[i].get(); }
}

// construct with root function
ODEData::ODEData(std::shared_ptr<ModPiece> rhs, std::shared_ptr<ModPiece> root, ref_vector<Eigen::VectorXd> const& refinputs, bool const autonomous, int const wrtIn) : rhs(rhs), root(root), autonomous(autonomous), wrtIn(wrtIn) {
  inputs.resize(refinputs.size());
  for( unsigned int i=0; i<refinputs.size(); ++i ) { inputs[i] = refinputs[i].get(); }
}

void ODEData::UpdateInputs(Eigen::VectorXd const& state, double const time) {
  if( autonomous ) {
    inputs[0] = state;
    return;
  }

  inputs[0] = Eigen::VectorXd::Constant(1, time);
  inputs[1] = state;
}
