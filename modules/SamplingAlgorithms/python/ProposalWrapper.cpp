#include "AllClassWrappers.h"

#include "AMProposal.h"
#include "InverseGammaProposal.h"
#include "MCMCProposal.h"
#include "MHProposal.h"

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
    
}