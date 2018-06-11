#include "AllClassWrappers.h"

#include "MUQ/SamplingAlgorithms/AMProposal.h"
#include "MUQ/SamplingAlgorithms/MCMCProposal.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::SamplingAlgorithms::PythonBindings;
namespace py = pybind11;

void muq::SamplingAlgorithms::PythonBindings::ProposalWrapper(py::module &m)
{
  py::class_<MCMCProposal, std::shared_ptr<MCMCProposal>> mcmcPro(m, "MCMCProposal");
  mcmcPro
    .def("Sample", &MCMCProposal::Sample)
    .def("LogDensity", &MCMCProposal::LogDensity);
    //.def("Construct", &MCMCProposal::Construct)
    //.def("GetMCMCProposalMap", &MCMCProposal::GetMCMCProposalMap);
  
  py::class_<MHProposal, MCMCProposal, std::shared_ptr<MHProposal>> mhPro(m, "MHProposal");
  mhPro
    .def(py::init<boost::property_tree::ptree const&, std::shared_ptr<AbstractSamplingProblem>>());
  
  py::class_<AMProposal, MHProposal, std::shared_ptr<AMProposal>> amPro(m, "AMProposal");
  amPro
    .def(py::init<boost::property_tree::ptree const&, std::shared_ptr<AbstractSamplingProblem>>())
    .def("Adapt", &AMProposal::Adapt);

}