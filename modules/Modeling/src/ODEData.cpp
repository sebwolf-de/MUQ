#include "MUQ/Modeling/ODEData.h"

#include "boost/none.hpp"

using namespace muq::Modeling;

// construct basic ode data
ODEData::ODEData(std::shared_ptr<ModPiece> rhs, ref_vector<Eigen::VectorXd> const& refinputs, bool const autonomous, int const wrtIn, int const wrtOut) : rhs(rhs), autonomous(autonomous), wrtIn(wrtIn), wrtOut(wrtOut) {
  inputs.resize(refinputs.size());
  for( unsigned int i=0; i<refinputs.size(); ++i ) { inputs[i] = refinputs[i].get(); }
}

// construct with root function
ODEData::ODEData(std::shared_ptr<ModPiece> rhs, std::shared_ptr<ModPiece> root, ref_vector<Eigen::VectorXd> const& refinputs, bool const autonomous, int const wrtIn, int const wrtOut) : rhs(rhs), root(root), autonomous(autonomous), wrtIn(wrtIn), wrtOut(wrtOut) {
  inputs.resize(refinputs.size());
  for( unsigned int i=0; i<refinputs.size(); ++i ) { inputs[i] = refinputs[i].get(); }
}
