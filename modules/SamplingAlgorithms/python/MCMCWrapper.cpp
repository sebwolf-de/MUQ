#include "AllClassWrappers.h"

#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;
namespace py = pybind11;

void PythonBindings::MCMCWrapper(py::module &m) {
  py::class_<SamplingAlgorithm, std::shared_ptr<SamplingAlgorithm>> sampAlg(m, "SamplingAlgorithm");
  sampAlg
    .def("Run", (SampleCollection const&  (SamplingAlgorithm::*)()) &SamplingAlgorithm::Run)
    .def("Run", (SampleCollection const&  (SamplingAlgorithm::*)(Eigen::VectorXd const&)) &SamplingAlgorithm::Run)
    .def("Run", (SampleCollection const&  (SamplingAlgorithm::*)(std::vector<Eigen::VectorXd> const&)) &SamplingAlgorithm::Run);

  py::class_<SingleChainMCMC, SamplingAlgorithm, std::shared_ptr<SingleChainMCMC>> singleMCMC(m, "SingleChainMCMC");
  singleMCMC
    .def("__init__", [](SingleChainMCMC &instance, py::dict d, std::shared_ptr<AbstractSamplingProblem> problem) {new (&instance) SingleChainMCMC(ConvertDictToPtree(d), problem);})
    .def("Kernels", &SingleChainMCMC::Kernels)
    .def("RunImpl", &SingleChainMCMC::RunImpl);
}
