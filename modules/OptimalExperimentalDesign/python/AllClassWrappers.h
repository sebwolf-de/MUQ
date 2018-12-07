#ifndef OPTIMALEXPERIMENTALDESIGN_ALLCLASSWRAPPERS_H_
#define OPTIMALEXPERIMENTALDESIGN_ALLCLASSWRAPPERS_H_

#include <pybind11/pybind11.h>

namespace muq {
  namespace OptimalExperimentalDesign {
    namespace PythonBindings {
      void EvidenceWrapper(pybind11::module &m);
      void OEDResidualWrapper(pybind11::module &m);
      void UtilityWrapper(pybind11::module &m);
    }
  }
}

#endif
