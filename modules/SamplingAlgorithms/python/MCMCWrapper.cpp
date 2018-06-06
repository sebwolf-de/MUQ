#include "AllClassWrappers.h"

#include "MCMC.h"
#include "MonteCarlo.h"
#include "SingleChainMCMC.h"

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
    
}