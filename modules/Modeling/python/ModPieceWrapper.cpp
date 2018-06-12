#include "AllClassWrappers.h"

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/IdentityOperator.h"
#include "MUQ/Modeling/ReplicateOperator.h"
#include "MUQ/Modeling/ConstantVector.h"
#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/ModGraphPiece.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

//using namespace muq::Modeling::PythonBindings;
using namespace muq::Modeling;
namespace py = pybind11;

class PyModPiece : public ModPiece {
public:
  /* Inherit the constructors */
  using ModPiece::ModPiece;

  /* Trampoline (need one for each virtual function) */
  void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) override {
    std::cout << "BLAH!" << std::endl;
    PYBIND11_OVERLOAD_PURE(
      void,          /* Return type */
      ModPiece,      /* Parent class */
      EvaluateImpl,  /* Name of function in C++ (must match Python name) */
      input          /* Argument(s) */
    );
    std::cout << "BLAH! BLAH!" << std::endl;
  }
};

class Publicist : public ModPiece {
public:
    using ModPiece::EvaluateImpl;
};

void muq::Modeling::PythonBindings::ModPieceWrapper(py::module &m)
{
    // Define some functions from the WorkPiece base class
    py::class_<ModPiece, PyModPiece, WorkPiece, std::shared_ptr<ModPiece>> mp(m, "ModPiece");
    mp
      .def(py::init<Eigen::VectorXi const&, Eigen::VectorXi const&>())
      .def("EvaluateImpl", (void (ModPiece::*)(ref_vector<Eigen::VectorXd> const&)) &Publicist::EvaluateImpl)
      .def("Evaluate", (std::vector<Eigen::VectorXd> const& (ModPiece::*)(std::vector<Eigen::VectorXd> const&)) &ModPiece::Evaluate)
      .def("Evaluate", (std::vector<Eigen::VectorXd> const& (ModPiece::*)()) &ModPiece::Evaluate)
      .def_readonly("inputSizes", &ModPiece::inputSizes)
      .def_readonly("outputSizes", &ModPiece::outputSizes)
      .def("GetRunTime", &ModPiece::GetRunTime)
      .def("ResetCallTime", &ModPiece::ResetCallTime)
      .def("GetNumCalls", &ModPiece::GetNumCalls)
      .def("Gradient", (Eigen::VectorXd const& (ModPiece::*)(unsigned int, unsigned int, std::vector<Eigen::VectorXd> const&, Eigen::VectorXd const&)) &ModPiece::Gradient)
      .def("Jacobian", (Eigen::MatrixXd const& (ModPiece::*)(unsigned int, unsigned int, std::vector<Eigen::VectorXd> const&)) &ModPiece::Jacobian)
      .def("ApplyJacobian", (Eigen::VectorXd const& (ModPiece::*)(unsigned int, unsigned int, std::vector<Eigen::VectorXd> const&, Eigen::VectorXd const&)) &ModPiece::ApplyJacobian)
      .def("GradientByFD", (Eigen::VectorXd (ModPiece::*)(unsigned int, unsigned int, std::vector<Eigen::VectorXd> const&, Eigen::VectorXd const&)) &ModPiece::GradientByFD)
      .def("JacobianByFD", (Eigen::MatrixXd (ModPiece::*)(unsigned int, unsigned int, std::vector<Eigen::VectorXd> const&)) &ModPiece::JacobianByFD)
      .def("ApplyJacobianByFD", (Eigen::VectorXd (ModPiece::*)(unsigned int, unsigned int, std::vector<Eigen::VectorXd> const&, Eigen::VectorXd const&)) &ModPiece::ApplyJacobianByFD);

    py::class_<ConstantVector, ModPiece, std::shared_ptr<ConstantVector>> cv(m, "ConstantVector");
    cv
        .def(py::init<Eigen::VectorXd const&>())
        .def("SetValue",  &ConstantVector::SetValue);

    py::class_<IdentityOperator, ModPiece, std::shared_ptr<IdentityOperator>> io(m, "IdentityOperator");
    io
        .def(py::init<unsigned int>());

    py::class_<ReplicateOperator, ModPiece, std::shared_ptr<ReplicateOperator>> ro(m, "ReplicateOperator");
    ro
      .def(py::init<unsigned int, unsigned int>());

    py::class_<ModGraphPiece, ModPiece, std::shared_ptr<ModGraphPiece>> mgp(m, "ModGraphPiece");
    mgp
      .def("GetGraph", &ModGraphPiece::GetGraph)
      .def("GetConstantPieces", &ModGraphPiece::GetConstantPieces);
}
