#ifndef PYDICTCONVERSION_H
#define PYDICTCONVERSION_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace muq{
  namespace Utilities{

    /** This function is useful for definining python wrappers of
    functions with ptree arguments.  It can be used in the python bindings
    to convert a python dictionary object into a ptree before calling the
    c++ function.  For example, it may be used like:
    @code
    struct PtreePrinter{
      void PrintPtree(boost::property_tree::ptree pt){
        boost::property_tree::write_json(std::cout, pt, true);
      }
    };

    // ... other stuff and module initialization

    py::class_<PtreePrinter, std::shared_ptr<PtreePrinter>> pp(m, "PtreePrinter");
    pp
      .def(py::init<>())
      .def("PrintPtree", [](PtreePrinter* p, py::dict d){p->PrintPtree(muq::Utilities::ConvertDictToPtree(d));});

    @endcode
    */
    boost::property_tree::ptree ConvertDictToPtree(pybind11::dict dict)
    {
      boost::property_tree::ptree pt;
      pybind11::object keys = pybind11::list(dict.attr("keys")());
      std::vector<std::string> keysCpp = keys.cast<std::vector<std::string>>();

      for(auto& key : keysCpp)
        pt.put(key,pybind11::str(dict.attr("get")(key)));

      return pt;
    }

  }
}


#endif
