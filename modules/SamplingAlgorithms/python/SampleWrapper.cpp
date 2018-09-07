#include "AllClassWrappers.h"

#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SampleCollection.h"
#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SamplingState.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::SamplingAlgorithms;
namespace py = pybind11;

void PythonBindings::SampleWrapper(py::module &m)
{
  py::class_<AbstractSamplingProblem, std::shared_ptr<AbstractSamplingProblem>> absSamp(m, "AbstractSamplingProblem");
  absSamp
    .def("LogDensity", &AbstractSamplingProblem::LogDensity)
    .def_readonly("numBlocks", &AbstractSamplingProblem::numBlocks)
    .def_readonly("blockSizes", &AbstractSamplingProblem::blockSizes);

  py::class_<SamplingProblem, AbstractSamplingProblem, std::shared_ptr<SamplingProblem>> sampProb(m, "SamplingProblem");
  sampProb
    .def(py::init<std::shared_ptr<muq::Modeling::ModPiece>>())
    .def("LogDensity", &SamplingProblem::LogDensity)
    .def("GradLogDensity", &SamplingProblem::GradLogDensity)
    .def("GetDistribution", &SamplingProblem::GetDistribution);

  py::class_<SamplingStateIdentity, std::shared_ptr<SamplingStateIdentity>> ssID(m, "SamplingStateIdentity");
  ssID
    .def(py::init<int>())
    .def_readonly("blockInd", &SamplingStateIdentity::blockInd);

  // py::class_<SamplingStatePartialMoment, std::shared_ptr<SamplingStatePartialMoment>> ssParMom(m, "SamplingStatePartialMoment");
  // ssParMom
  //   .def(py::init<int, int, Eigen::VectorXd const&>())
  //   .def_readonly("blockInd", &SamplingStatePartialMoment::blockInd)
  //   .def_readonly("momentPower", &SamplingStatePartialMoment::momentPower);
  //   //.def_readonly("mu", &SamplingStatePartialMoment::mu);

  py::class_<SampleCollection, std::shared_ptr<SampleCollection>> sampColl(m, "SampleCollection");
  sampColl
    .def("__getitem__", (const std::shared_ptr<SamplingState> (SampleCollection::*)(unsigned) const) &SampleCollection::at)
    .def("size", &SampleCollection::size)
    .def("CentralMoment", (Eigen::VectorXd (SampleCollection::*)(unsigned, int) const) &SampleCollection::CentralMoment, py::arg("order"), py::arg("blockDim") = -1)
    .def("CentralMoment", (Eigen::VectorXd (SampleCollection::*)(unsigned, Eigen::VectorXd const&, int) const) &SampleCollection::CentralMoment, py::arg("order"), py::arg("mean"), py::arg("blockDim") = -1)
    .def("Mean", &SampleCollection::Mean, py::arg("blockDim") = -1)
    .def("Variance", &SampleCollection::Variance, py::arg("blockDim") = -1)
    .def("Covariance", (Eigen::MatrixXd (SampleCollection::*)(int) const) &SampleCollection::Covariance, py::arg("blockDim") = -1)
    .def("Covariance", (Eigen::MatrixXd (SampleCollection::*)(Eigen::VectorXd const&, int) const) &SampleCollection::Covariance, py::arg("mean"), py::arg("blockDim") = -1)
    .def("ESS", &SampleCollection::ESS, py::arg("blockDim")=-1)
    .def("Add", &SampleCollection::Add)
    .def("Weights", &SampleCollection::Weights)
    .def("AsMatrix", &SampleCollection::AsMatrix, py::arg("blockDim")=-1);

  py::class_<SamplingState, std::shared_ptr<SamplingState>> sampState(m, "SamplingState");
  sampState
    .def(py::init<Eigen::VectorXd const&>())
    .def(py::init<Eigen::VectorXd const&, double>())
    .def(py::init<std::vector<Eigen::VectorXd> const&>())
    .def(py::init<std::vector<Eigen::VectorXd> const&, double>())
    .def_readonly("weight", &SamplingState::weight)
    .def_readonly("state", &SamplingState::state)
    .def("HasMeta", &SamplingState::HasMeta)
    .def("TotalDim", &SamplingState::TotalDim);

}
