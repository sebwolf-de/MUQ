#include "MUQ/Modeling/ODEData.h"

#include "boost/none.hpp"

using namespace muq::Modeling;

ODEData::ODEData(std::shared_ptr<WorkPiece> rhs, ref_vector<boost::any> const& inputs) : rhs(rhs), inputs(inputs) {}

ODEData::ODEData(std::shared_ptr<WorkPiece> rhs, std::shared_ptr<WorkPiece> root, ref_vector<boost::any> const& inputs) : rhs(rhs), root(root), inputs(inputs) {}

