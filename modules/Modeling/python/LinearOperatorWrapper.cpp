#include "AllClassWrappers.h"

#include "MUQ/Modeling/ModPiece.h"

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/EigenLinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Modeling::PythonBindings;
using namespace muq::Modeling;
namespace py = pybind11;


void muq::Modeling::PythonBindings::LinearOperatorWrapper(py::module &m)
{

  py::class_<LinearOperator, ModPiece, WorkPiece, std::shared_ptr<LinearOperator>> lo(m, "LinearOperator");
  lo
    .def("rows", &LinearOperator::rows)
    .def("cols", &LinearOperator::cols)
    .def("GetMatrix", &LinearOperator::GetMatrix);

  py::class_<EigenLinearOperator<Eigen::MatrixXd>, LinearOperator, ModPiece, WorkPiece, std::shared_ptr<EigenLinearOperator<Eigen::MatrixXd>>> elo(m, "DenseLinearOperator");
  elo
    .def(py::init<Eigen::MatrixXd>())
    .def("Apply", &EigenLinearOperator<Eigen::MatrixXd>::Apply)
    .def("ApplyTranspose", &EigenLinearOperator<Eigen::MatrixXd>::ApplyTranspose);

  py::class_<IdentityOperator, LinearOperator, ModPiece, WorkPiece, std::shared_ptr<IdentityOperator>> io(m, "IdentityOperator");
  io
    .def(py::init<unsigned int>())
    .def("Apply", &IdentityOperator::Apply)
    .def("ApplyTranspose", &IdentityOperator::ApplyTranspose);

};
