#include "AllClassWrappers.h"

//#include "MUQ/SamplingAlgorithms/ISKernel.h"
//#include "MUQ/SamplingAlgorithms/MCKernel.h"
#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::SamplingAlgorithms::PythonBindings;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;
namespace py = pybind11;

std::shared_ptr<MHKernel> CreatePythonMHKernel(py::dict d, std::shared_ptr<AbstractSamplingProblem> problem)
{
    return std::make_shared<MHKernel>(ConvertDictToPtree(d), problem);
};


void muq::SamplingAlgorithms::PythonBindings::KernelWrapper(py::module &m)
{
  
  
  py::class_<TransitionKernel, std::shared_ptr<TransitionKernel>> transKern(m, "TransitionKernel");
  transKern
  //  .def("Construct", &TransitionKernel::Construct)
  //  .def("GetTransitionKernelMap", &TransitionKernel::GetTransitionKernelMap)
    .def("PreStep", &TransitionKernel::PreStep)
    .def("PostStep", &TransitionKernel::PostStep)
    .def("Step", &TransitionKernel::Step)
    .def_readonly("blockInd", &TransitionKernel::blockInd);

   
  py::class_<MHKernel, TransitionKernel, std::shared_ptr<MHKernel>> mhKern(m, "MHKernel");
  mhKern
    .def(py::init<>(&CreatePythonMHKernel))
    .def("Proposal", &MHKernel::Proposal)
    .def("PostStep", &MHKernel::PostStep)
    .def("Step", &MHKernel::Step);

}