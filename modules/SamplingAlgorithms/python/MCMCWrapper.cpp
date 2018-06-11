#include "AllClassWrappers.h"

//#include "MUQ/SamplingAlgorithms/MonteCarlo.h"
#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::SamplingAlgorithms::PythonBindings;
namespace py = pybind11;

void muq::SamplingAlgorithms::PythonBindings::MCMCWrapper(py::module &m)
{
  py::class_<SamplingAlgorithm, std::shared_ptr<SamplingAlgorithm>> sampAlg(m, "SamplingAlgorithm");
  sampAlg
    .def("Run", (SampleCollection const&  (SamplingAlgorithm::*)()) &SamplingAlgorithm::Run)
    .def("Run", (SampleCollection const&  (SamplingAlgorithm::*)(Eigen::VectorXd const&)) &SamplingAlgorithm::Run)
    .def("Run", (SampleCollection const&  (SamplingAlgorithm::*)(std::vector<Eigen::VectorXd> const&)) &SamplingAlgorithm::Run);

  py::class_<SingleChainMCMC, SamplingAlgorithm, std::shared_ptr<SingleChainMCMC>> singleMCMC(m, "SingleChainMCMC");
  singleMCMC
    .def("Kernels", &SingleChainMCMC::Kernels)
    .def("RunImpl", &SingleChainMCMC::RunImpl);
}
