#include "MUQ/Modeling/ModPiece.h"

#include <pybind11/pybind11.h>

using namespace muq::Modeling;
namespace py = pybind11;

class PyModPiece : public ModPiece {
public:
  /* Inherit the constructors */
  using ModPiece::ModPiece;

  /* Trampoline (need one for each virtual function) */
  void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) override {
    PYBIND11_OVERLOAD_PURE(
      void,          /* Return type */
      ModPiece,      /* Parent class */
      EvaluateImpl,  /* Name of function in C++ (must match Python name) */
      input          /* Argument(s) */
    );
  }
};
