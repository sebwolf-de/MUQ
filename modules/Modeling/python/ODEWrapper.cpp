#include "AllClassWrappers.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include "MUQ/Modeling/ODEPiece.h"

using namespace muq::Modeling;
using namespace muq::Utilities;
namespace py = pybind11;

void PythonBindings::ODEWrapper(py::module &m) {
  py::class_<ODEPiece, ModPiece, std::shared_ptr<ODEPiece> > ode(m, "ODEPiece");
  ode.def(py::init( [] (std::shared_ptr<ModPiece> const& rhs, py::dict const& d) { return new ODEPiece(rhs, ConvertDictToPtree(d)); }), py::keep_alive<1, 2>());

}
