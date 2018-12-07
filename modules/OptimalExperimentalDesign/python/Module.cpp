#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "AllClassWrappers.h"

using namespace muq::OptimalExperimentalDesign::PythonBindings;
namespace py = pybind11;

PYBIND11_MODULE(pymuqOptimalExperimentalDesign, m) {
  EvidenceWrapper(m);
  OEDResidualWrapper(m);
  UtilityWrapper(m);
}
