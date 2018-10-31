#include "AllClassWrappers.h"

#include "MUQ/Approximation/Quadrature/GaussQuadrature.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Approximation::PythonBindings;
using namespace muq::Modeling;

namespace py = pybind11;

void muq::Approximation::PythonBindings::QuadratureWrapper(py::module &m)
{

  py::class_<GaussQuadrature, std::shared_ptr<GaussQuadrature>> quadFunc(m, "GaussQuadrature");
  quadFunc
    .def(py::init<std::shared_ptr<OrthogonalPolynomial>, int>())
    .def("Compute", &GaussQuadrature::Compute)
    .def("Points", &GaussQuadrature::Points)
    .def("Weights", &GaussQuadrature::Weights);

}
