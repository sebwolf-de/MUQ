#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Modeling/IdentityPiece.h"
#include "MUQ/Modeling/ConstantPiece.h"

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/eigen.h"

#include <string>

#include <boost/any.hpp>
#include "PyAny.h"

using namespace muq::Modeling;
namespace py = pybind11;


void AnyFunction(boost::any x)
{
    std::cout << "Success!  Ran AnyFunction with type " << x.type().name() << std::endl;
}


PYBIND11_PLUGIN(muqModeling) {
    py::module m("muqModeling", "Python bindings for the muqModeling library.");

    // py::class_<boost::any> any(m, "BoostAny");
    // any
    //     .def(py::init<double>())
    //     .def(py::init<int>())
    //     .def(py::init<std::string>())
    //     .def(py::init<Eigen::Ref<Eigen::MatrixXd>>())
    //     .def("Get", 
    
    // py::implicitly_convertible<int, boost::any>();
    // py::implicitly_convertible<double, boost::any>();
    // py::implicitly_convertible<std::string, boost::any>();
    // py::implicitly_convertible<Eigen::Ref<Eigen::MatrixXd>, boost::any>();
    
    m.def("AnyFunction", &AnyFunction);

    // Define some functions from the WorkPiece base class
    py::class_<WorkPiece> wp(m, "WorkPiece");
    wp
        .def("Evaluate", (std::vector<boost::any> (WorkPiece::*)(std::vector<boost::any> const&)) &WorkPiece::Evaluate)
        .def("Evaluate", (std::vector<boost::any> (WorkPiece::*)()) &WorkPiece::Evaluate);
         
    // The Identity WorkPiece
    py::class_<IdentityPiece> ip(m, "IdentityPiece", wp);
    ip
        .def(py::init())
        .def(py::init<int const>());

    // The Constant WorkPiece
    py::class_<ConstantPiece> cp(m, "ConstantPiece", wp);
    cp
        .def(py::init<std::vector<boost::any> const&>())
        .def(py::init<boost::any const&>());
    
    
    return m.ptr();
}
