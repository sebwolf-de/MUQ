#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Modeling/IdentityPiece.h"
#include "MUQ/Modeling/ConstantPiece.h"
#include "MUQ/Modeling/WorkGraphPiece.h"
#include "MUQ/Modeling/WorkGraph.h"

#include "MUQ/Modeling/Dolfin/FenicsPiece.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <boost/any.hpp>
#include <functional>
#include <vector>

#include "PyAny.h"

using namespace muq::Modeling;
namespace py = pybind11;


void AnyFunction(boost::any x)
{
    std::cout << "Success!  Ran AnyFunction with type " << x.type().name() << std::endl;
}

PYBIND11_PLUGIN(pymuqModeling) {
    py::module m("pymuqModeling", "Python bindings for the muqModeling library.");

    m.def("AnyFunction", &AnyFunction);
    
    // Define some functions from the WorkPiece base class
    py::class_<WorkPiece,  std::shared_ptr<WorkPiece>> wp(m, "WorkPiece");
    wp
        .def("Evaluate", (std::vector<boost::any> (WorkPiece::*)(std::vector<boost::any> const&)) &WorkPiece::Evaluate)
        .def("Evaluate", (std::vector<boost::any> (WorkPiece::*)()) &WorkPiece::Evaluate);
         
    // The Identity WorkPiece
    py::class_<IdentityPiece, std::shared_ptr<IdentityPiece>> ip(m, "IdentityPiece", wp);
    ip
        .def(py::init())
        .def(py::init<int const>())
        .def(py::init<std::vector<std::string> const&>())
        .def(py::init<std::map<unsigned int, std::string> const&>())
        .def(py::init<std::map<unsigned int, std::string> const&, unsigned int const>());

    // The Constant WorkPiece
    py::class_<ConstantPiece, std::shared_ptr<ConstantPiece>> cp(m, "ConstantPiece", wp);
    cp
        .def(py::init<std::vector<boost::any> const&>())
        .def(py::init<boost::any const&>());
    
    // The WorkGraphPiece
    py::class_<WorkGraphPiece, std::shared_ptr<WorkGraphPiece>> wgp(m, "WorkGraphPiece", wp);
           
    // The WorkGraph
    py::class_<WorkGraph> wg(m, "WorkGraph");
    wg
        .def(py::init())
        .def("NumNodes", &WorkGraph::NumNodes)
        .def("NumEdges", &WorkGraph::NumEdges)
        .def("AddNode", &WorkGraph::AddNode)
        .def("AddEdge", &WorkGraph::AddEdge)
        .def("HasNode", (bool (WorkGraph::*)(std::string const&) const) &WorkGraph::HasNode)
        .def("Visualize", &WorkGraph::Visualize)
        .def("DependentCut", &WorkGraph::DependentCut)
        .def("CreateWorkPiece", &WorkGraph::CreateWorkPiece)
        .def("Constant", (bool (WorkGraph::*)(std::string const&) const) &WorkGraph::Constant)
        .def("GetConstantOutputs", (std::vector<boost::any> const& (WorkGraph::*)(std::string const&) const) &WorkGraph::GetConstantOutputs);
    

    py::class_<FenicsPiece, std::shared_ptr<FenicsPiece>> fp(m, "FenicsPiece", wp);
    fp
        .def(py::init<pybind11::object const&, pybind11::object const&, std::vector<pybind11::object> const&>() )
        .def("EvaluateVec", &FenicsPiece::EvaluateVec);
    
    return m.ptr();
}
