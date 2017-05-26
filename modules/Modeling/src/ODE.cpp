#include "MUQ/Modeling/ODE.h"

using namespace muq::Modeling;

ODE::ODE(std::shared_ptr<WorkPiece> rhs) : ODEBase(rhs) {}

ODE::~ODE() {}

void ODE::EvaluateImpl(ref_vector<boost::any> const& inputs) {
}
