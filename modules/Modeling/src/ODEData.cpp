#include "MUQ/Modeling/ODEData.h"

#include "boost/none.hpp"

using namespace muq::Modeling;

// construct basic ode data
ODEData::ODEData(std::shared_ptr<WorkPiece> rhs, ref_vector<boost::any> const& inputs, bool const autonomous, int const wrtIn, int const wrtOut) : rhs(rhs), inputs(inputs), autonomous(autonomous), wrtIn(wrtIn), wrtOut(wrtOut) {}

// construct with root function
ODEData::ODEData(std::shared_ptr<WorkPiece> rhs, std::shared_ptr<WorkPiece> root, ref_vector<boost::any> const& inputs, bool const autonomous, int const wrtIn, int const wrtOut) : rhs(rhs), root(root), inputs(inputs), autonomous(autonomous), wrtIn(wrtIn), wrtOut(wrtOut) {}

