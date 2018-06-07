#include "AllClassWrappers.h"

#include "MUQ/Utilities/MultiIndices/MultiIndex.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexLimiter.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexSet.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::SamplingAlgorithms::PythonBindings;
namespace py = pybind11;

void muq::SamplingAlgorithms::PythonBindings::MultiIndicesWrapper(py::module &m)
{
    
}