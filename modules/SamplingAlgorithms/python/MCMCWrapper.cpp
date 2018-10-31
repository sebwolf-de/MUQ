#include "AllClassWrappers.h"

#include "MUQ/SamplingAlgorithms/SingleChainMCMC.h"
#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"
#include "MUQ/SamplingAlgorithms/MCMCFactory.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;
namespace py = pybind11;

void PythonBindings::MCMCWrapper(py::module &m) {
  py::class_<SamplingAlgorithm, std::shared_ptr<SamplingAlgorithm>> sampAlg(m, "SamplingAlgorithm");
  sampAlg
    .def("Run", (std::shared_ptr<SampleCollection>  (SamplingAlgorithm::*)()) &SamplingAlgorithm::Run,
                 py::call_guard<py::scoped_ostream_redirect,py::scoped_estream_redirect>())
    .def("Run", (std::shared_ptr<SampleCollection>  (SamplingAlgorithm::*)(Eigen::VectorXd const&)) &SamplingAlgorithm::Run,
                 py::call_guard<py::scoped_ostream_redirect,py::scoped_estream_redirect>())
    .def("Run", (std::shared_ptr<SampleCollection>  (SamplingAlgorithm::*)(std::vector<Eigen::VectorXd> const&)) &SamplingAlgorithm::Run,
                 py::call_guard<py::scoped_ostream_redirect,py::scoped_estream_redirect>())
    .def("GetSamples", &SamplingAlgorithm::GetSamples);

  py::class_<SingleChainMCMC, SamplingAlgorithm, std::shared_ptr<SingleChainMCMC>> singleMCMC(m, "SingleChainMCMC");
  singleMCMC
    .def(py::init( [](py::dict d, std::shared_ptr<AbstractSamplingProblem> problem) {return new SingleChainMCMC(ConvertDictToPtree(d), problem);}))
    .def(py::init( [](py::dict d, std::shared_ptr<AbstractSamplingProblem> problem, std::vector<std::shared_ptr<TransitionKernel>> kernelsIn) {return new SingleChainMCMC(ConvertDictToPtree(d), problem, kernelsIn);}))
    .def("Kernels", &SingleChainMCMC::Kernels)
    .def("RunImpl", &SingleChainMCMC::RunImpl);

  py::class_<MCMCFactory, std::shared_ptr<MCMCFactory>> fact(m, "MCMCFactory");
  fact
    .def_static("CreateSingleChain", &MCMCFactory::CreateSingleChain);
}
