#include "AllClassWrappers.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include "MUQ/OptimalExperimentalDesign/Utility.h"

namespace py = pybind11;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::OptimalExperimentalDesign;

void PythonBindings::UtilityWrapper(py::module &m) {
  py::class_<Utility, ModPiece, std::shared_ptr<Utility> > util(m, "Utility");

  util.def(py::init( [] (std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, py::dict d) { return new Utility(prior, likelihood, evidence, ConvertDictToPtree(d)); }));

  util.def(py::init( [] (std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, std::shared_ptr<Distribution> const& biasing, py::dict d) { return new Utility(prior, likelihood, evidence, biasing, ConvertDictToPtree(d)); }));

  util.def(py::init( [] (std::shared_ptr<Distribution> const& prior, std::shared_ptr<OEDResidual> const& resid, py::dict d) { return new Utility(prior, resid, ConvertDictToPtree(d)); }));

  util.def(py::init( [] (std::shared_ptr<Distribution> const& prior, std::shared_ptr<OEDResidual> const& resid, std::shared_ptr<Distribution> const& biasing, py::dict d) { return new Utility(prior, resid, biasing, ConvertDictToPtree(d)); }));

#if MUQ_HAS_PARCER==1
  util.def(py::init( [] (std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, py::dict d, std::shared_ptr<parcer::Communicator> const& comm) { return new Utility(prior, likelihood, evidence, ConvertDictToPtree(d), comm); }));

  util.def(py::init( [] (std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, std::shared_ptr<Distribution> const& biasing, py::dict d, std::shared_ptr<parcer::Communicator> const& comm) { return new Utility(prior, likelihood, evidence, biasing, ConvertDictToPtree(d), comm); }));

  util.def(py::init( [] (std::shared_ptr<Distribution> const& prior, std::shared_ptr<OEDResidual> const& resid, py::dict d, std::shared_ptr<parcer::Communicator> const& comm) { return new Utility(prior, resid, ConvertDictToPtree(d), comm); }));

  util.def(py::init( [] (std::shared_ptr<Distribution> const& prior, std::shared_ptr<OEDResidual> const& resid, std::shared_ptr<Distribution> const& biasing, py::dict d, std::shared_ptr<parcer::Communicator> const& comm) { return new Utility(prior, resid, biasing, ConvertDictToPtree(d), comm); }));
#endif
}
