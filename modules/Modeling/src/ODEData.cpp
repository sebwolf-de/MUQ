#include "MUQ/Modeling/ODEData.h"

#include "boost/none.hpp"

using namespace muq::Modeling;

// construct basic ode data
ODEData::ODEData(std::shared_ptr<WorkPiece> rhs, ref_vector<boost::any> const& inputs, int const wrtIn, int const wrtOut) : rhs(rhs), inputs(inputs), wrtIn(wrtIn), wrtOut(wrtOut) {}

// construct with root function
ODEData::ODEData(std::shared_ptr<WorkPiece> rhs, std::shared_ptr<WorkPiece> root, ref_vector<boost::any> const& inputs, int const wrtIn, int const wrtOut) : rhs(rhs), root(root), inputs(inputs), wrtIn(wrtIn), wrtOut(wrtOut) {}

