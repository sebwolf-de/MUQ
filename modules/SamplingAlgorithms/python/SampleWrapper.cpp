#include "AllClassWrappers.h"

#include "AbstractSamplingProblem.h"
#include "ImportanceSampling.h"
#include "SampleCollection.h"
#include "SamplingAlgorithm.h"
#include "SamplingProblem.h"
#include "SamplingState.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::SamplingAlgorithms::PythonBindings;
namespace py = pybind11;

void muq::SamplingAlgorithms::PythonBindings::SampleWrapper(py::module &m)
{
  py::class_<AbstractSamplingProblem, std::shared_ptr<AbstractSamplingProblem>> absSamp(m, "AbstractSamplingProblem");
  absSamp
    .def("LogDensity", &AbstractSamplingProblem::LogDensity)
    .def("GradLogDensity", &AbstractSamplingProblem::GradLogDensity)
    .def_readonly("numBlocks", &AbstractSamplingProblem::numBlocks)
    .def_readonly("blockSizes", &AbstractSamplingProblem::blockSizes);
  /*  
  py::class_<SamplingAlgorithm, std::shared_ptr<SamplingAlgorithm>> sampAlg(m, "SamplingAlgorithm");
  sampAlg 
    .def("GetSamples", &SamplingAlgorithms::GetSamples)
    .def("Run", (SampleCollection const& (SamplingAlgorithm::*)()) &SamplingAlgorithm::Run)
    .def("Run", (SampleCollection const& (SamplingAlgorithm::*)(Eigen::VectorXd const&)) &SamplingAlgorithm::Run)
    .def("Run", (SampleCollection const& (SamplingAlgorithm::*)(std::vector<Eigen::VectorXd> const& )) &SamplingAlgorithm::Run)
    .def("RunImpl", &SamplingAlgorithms::RunImpl);
    */
    
  py::class_<SamplingProblem, std::shared_ptr<SamplingProblem>> sampProb(m, "SamplingProblem");
  SamplingProblem
    .def("LogDensity", &SamplingProblems::LogDensity)
    .def("GradLogDensity", &SamplingProblems::GradLogDensity)
    .def("GetNumBlocks", &SamplingProblems::GeNumBlocks)
    .def("GetBlockSizes", &SamplingProblems::GetBlockSizes);
}