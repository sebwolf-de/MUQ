#include "AllClassWrappers.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include "MUQ/OptimalExperimentalDesign/Evidence.h"

namespace py = pybind11;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::OptimalExperimentalDesign;

void PythonBindings::EvidenceWrapper(py::module &m) {
  py::class_<Evidence, Distribution, std::shared_ptr<Evidence> > evi(m, "Evidence");

  evi.def(py::init( [] (std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, py::dict const& d) { return new Evidence(prior, likelihood, ConvertDictToPtree(d)); }));

  evi.def(py::init( [] (std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& biasing, py::dict const& d) { return new Evidence(prior, likelihood, biasing, ConvertDictToPtree(d)); }));

#if MUQ_HAS_PARCER==1
  evi.def(py::init( [] (std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, py::dict const& d, std::shared_ptr<parcer::Communicator> const& comm) { return new Evidence(prior, likelihood, ConvertDictToPtree(d), comm); }));

  evi.def(py::init( [] (std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& biasing, py::dict const& d, std::shared_ptr<parcer::Communicator> const& comm) { return new Evidence(prior, likelihood, biasing, ConvertDictToPtree(d), comm); }));
#endif
}
