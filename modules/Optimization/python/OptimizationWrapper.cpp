#include "AllClassWrappers.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>

#include <string>

#include <functional>
#include <vector>

#include "MUQ/Optimization/Optimization.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include "MUQ/Modeling/Python/PyAny.h"

using namespace muq::Utilities;
using namespace muq::Optimization;
namespace py = pybind11;

void PythonBindings::OptimizationWrapper(pybind11::module &m) {
  py::class_<Optimization, std::shared_ptr<Optimization> > opt(m, "Optimization");
  opt
  .def(py::init( [](std::shared_ptr<CostFunction> cost, py::dict d) { return new Optimization(cost, ConvertDictToPtree(d)); }))
  .def("AddInequalityConstraint", &Optimization::AddInequalityConstraint)
  .def("Solve", (std::pair<Eigen::VectorXd, double> (Optimization::*)(std::vector<Eigen::VectorXd> const&)) &Optimization::Solve);
}
