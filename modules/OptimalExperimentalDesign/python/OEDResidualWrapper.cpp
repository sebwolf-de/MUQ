#include "AllClassWrappers.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include "MUQ/OptimalExperimentalDesign/OEDResidual.h"

namespace py = pybind11;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::OptimalExperimentalDesign;

void PythonBindings::OEDResidualWrapper(py::module &m) {
  py::class_<OEDResidual, ModPiece, std::shared_ptr<OEDResidual> > res(m, "OEDResidual");

  res.def(py::init( [] (std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, py::dict d) { return new OEDResidual(likelihood, evidence, ConvertDictToPtree(d)); }));

#if MUQ_HAS_PARCER==1
  res.def(py::init( [] (std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, py::dict d, std::shared_ptr<parcer::Communicator> const& comm) { return new OEDResidual(likelihood, evidence, ConvertDictToPtree(d), comm); }));
#endif
}
