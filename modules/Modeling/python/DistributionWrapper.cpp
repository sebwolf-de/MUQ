#include "AllClassWrappers.h"

#include "MUQ/Modeling/ModPiece.h"

#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/Distributions/RandomVariable.h"
#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/InverseGamma.h"
#include "MUQ/Modeling/Distributions/UniformBox.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Modeling::PythonBindings;
namespace py = pybind11;


void muq::Modeling::PythonBindings::DistributionWrapper(py::module &m)
{
    py::class_<Distribution, std::shared_ptr<Distribution>> dist(m, "Distribution");
    dist
      .def("LogDensity", (double (Distribution::*)(std::vector<Eigen::VectorXd> const&)) &Distribution::LogDensity)
      .def("LogDensity", (double (Distribution::*)(Eigen::VectorXd const&)) &Distribution::LogDensity)
      .def("GradLogDensity", (Eigen::VectorXd (Distribution::*)(unsigned int, std::vector<Eigen::VectorXd> const&)) &Distribution::GradLogDensity)
      //.def("GradLogDensity", (Eigen::VectorXd (Distribution::*)(unsigned int, Eigen::VectorXd const&)) &Distribution::GradLogDensity)
      .def("Sample", (Eigen::VectorXd (Distribution::*)(std::vector<Eigen::VectorXd> const&)) &Distribution::Sample)
      //.def("Sample", (Eigen::VectorXd (Distribution::*)(Eigen::VectorXd const&)) &Distribution::Sample)
      .def("Sample", (Eigen::VectorXd (Distribution::*)()) &Distribution::Sample)
      .def("AsDensity", &Distribution::AsDensity)
      .def("AsVariable", &Distribution::AsVariable)
      .def_readonly("varSize", &Distribution::varSize)
      .def_readonly("hyperSizes", &Distribution::hyperSizes);

    py::class_<DensityBase, Distribution, ModPiece, std::shared_ptr<DensityBase>> densBase(m, "DensityBase");
    densBase
      .def(py::init<Eigen::VectorXi>());

    py::class_<Density, DensityBase, std::shared_ptr<Density>> dens(m, "Density");
    dens
      .def("GetDistribution", &Density::GetDistribution);

    py::class_<RandomVariable, Distribution, ModPiece, std::shared_ptr<RandomVariable>> rv(m, "RandomVariable");
    rv
      .def(py::init<std::shared_ptr<Distribution>>());

    py::class_<UniformBox, Distribution, std::shared_ptr<UniformBox>> uBox(m, "UniformBox");
    uBox
      .def(py::init<Eigen::MatrixXd const&>());

    py::class_<Gaussian, Distribution, std::shared_ptr<Gaussian>> gauss(m,"Gaussian");
    gauss
      .def(py::init<unsigned int>())
      .def(py::init<unsigned int, Gaussian::InputMask>())
      .def(py::init<Eigen::VectorXd const&>())
      .def(py::init<Eigen::VectorXd const&, Gaussian::InputMask>())
      .def(py::init<Eigen::VectorXd const&, Eigen::MatrixXd const&>())
      .def(py::init<Eigen::VectorXd const&, Eigen::MatrixXd const&, Gaussian::Mode>())
      .def(py::init<Eigen::VectorXd const&, Eigen::MatrixXd const&, Gaussian::Mode, Gaussian::InputMask>())
      .def("GetMode", &Gaussian::GetMode)
      .def("Dimension", &Gaussian::Dimension)
      .def("GetCovariance", &Gaussian::GetCovariance)
      .def("GetPrecision", &Gaussian::GetPrecision)
      .def("GetMean", &Gaussian::GetMean)
      .def("SetMean", &Gaussian::SetMean)
      .def("SetCovariance", &Gaussian::SetCovariance)
      .def("SetPrecision", &Gaussian::SetPrecision);

    py::enum_<Gaussian::Mode>(gauss, "Mode")
          .value("Covariance", Gaussian::Mode::Covariance)
          .value("Precision", Gaussian::Mode::Precision)
          .export_values();

    py::enum_<Gaussian::ExtraInputs>(gauss, "ExtraInputs", py::arithmetic())
          .value("None", Gaussian::ExtraInputs::None)
          .value("Mean", Gaussian::ExtraInputs::Mean)
          .value("DiagCovariance", Gaussian::ExtraInputs::DiagCovariance)
          .value("DiagPrecision", Gaussian::ExtraInputs::DiagPrecision)
          .value("FullCovariance", Gaussian::ExtraInputs::FullCovariance)
          .value("FullPrecision", Gaussian::ExtraInputs::FullPrecision)
          .export_values();


    py::class_<InverseGamma, Distribution, std::shared_ptr<InverseGamma>> ig(m, "InverseGamma");
    ig
      .def(py::init<double,double>())
      .def_readonly("alpha", &InverseGamma::alpha)
      .def_readonly("alpha", &InverseGamma::beta);
}
