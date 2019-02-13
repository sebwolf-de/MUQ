#include "AllClassWrappers.h"

#include "MUQ/Approximation/Quadrature/GaussQuadrature.h"
#include "MUQ/Approximation/Quadrature/ClenshawCurtisQuadrature.h"
#include "MUQ/Approximation/Quadrature/FullTensorQuadrature.h"
#include "MUQ/Approximation/Quadrature/SmolyakQuadrature.h"

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

  py::class_<Quadrature, std::shared_ptr<Quadrature>> quadBase(m, "Quadrature");
  quadBase
    .def("Compute",(void (Quadrature::*)(unsigned int)) &Quadrature::Compute)
    .def("Compute",(void (Quadrature::*)(Eigen::RowVectorXi const&)) &Quadrature::Compute)
    .def("Dim", &Quadrature::Dim)
    .def("Points", &Quadrature::Points)
    .def("Weights", &Quadrature::Weights);

  py::class_<GaussQuadrature, Quadrature, std::shared_ptr<GaussQuadrature>> quadFunc(m, "GaussQuadrature");
  quadFunc
    .def(py::init<std::shared_ptr<OrthogonalPolynomial>>())
    .def(py::init<std::shared_ptr<OrthogonalPolynomial>, int>())
    .def("Compute", &GaussQuadrature::Compute);

  py::class_<ClenshawCurtisQuadrature, Quadrature, std::shared_ptr<ClenshawCurtisQuadrature>> ccQuad(m,"ClenshawCurtisQuadrature");
  ccQuad
    .def(py::init<>())
    .def(py::init<bool>())
    .def("Compute", &ClenshawCurtisQuadrature::Compute);

  py::class_<FullTensorQuadrature, Quadrature, std::shared_ptr<FullTensorQuadrature>> tensQuad(m,"FullTensorQuadrature");
  tensQuad
    .def(py::init<unsigned int, std::shared_ptr<Quadrature>>())
    .def(py::init<unsigned int, std::shared_ptr<Quadrature>, unsigned int>())
    .def(py::init<std::vector<std::shared_ptr<Quadrature>>, Eigen::RowVectorXi>());

  py::class_<SmolyakQuadrature, Quadrature, std::shared_ptr<SmolyakQuadrature>> smolyQuad(m,"SmolyakQuadrature");
  smolyQuad
    .def(py::init<unsigned int, std::shared_ptr<Quadrature> const&>())
    .def(py::init<std::vector<std::shared_ptr<Quadrature>> const&>())
    .def("Compute", (void (SmolyakQuadrature::*)(unsigned int)) &SmolyakQuadrature::Compute)
    .def("Compute", (void (SmolyakQuadrature::*)(Eigen::RowVectorXi const&)) &SmolyakQuadrature::Compute)
    .def("Compute", (void (SmolyakQuadrature::*)(std::shared_ptr<muq::Utilities::MultiIndexSet> const&)) &SmolyakQuadrature::Compute)
    .def("ComputeWeights", &SmolyakQuadrature::ComputeWeights)
    .def("BuildMultis", &SmolyakQuadrature::BuildMultis);
}
