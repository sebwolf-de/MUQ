#include "AllClassWrappers.h"

#include "Attributes.h"
#include "BlockDataset.h"
#include "H5Object.h"
#include "HDF5File.h"
#include "HDF5Types.h"
#include "Pathtools.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::SamplingAlgorithms::PythonBindings;
namespace py = pybind11;

void muq::SamplingAlgorithms::PythonBindings::HDF5Wrapper(py::module &m)
{
    
}