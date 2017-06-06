#ifndef DOLFINPIECE_H_
#define DOLFINPIECE_H_

#include "MUQ/Modeling/WorkPiece.h"

#include <pybind11/pybind11.h>

namespace muq{
namespace Modeling{

    class DolfinPiece : public WorkPiece{
    
    public:

        DolfinPiece(pybind11::object const& problemIn, pybind11::object const& outputFieldIn, std::vector<pybind11::object> inputs);

        void EvaluateFunc(std::vector<pybind11::object> const& inputs);
    
        Eigen::VectorXd EvaluateVec(Eigen::Ref<const Eigen::VectorXd> const& x, std::vector<double> vals);
    
    private:

        static void CheckType(pybind11::object const& obj, std::string const& requiredType);

        void ExtractInputs(std::vector<pybind11::object> const& inputs);
    
        std::map<unsigned, std::pair<pybind11::object, std::vector<std::string>>> inputExprs;
        std::map<unsigned, std::shared_ptr<dolfin::Function>> inputFuncs;
    
        std::shared_ptr<dolfin::Function> outputField;
        std::shared_ptr<dolfin::LinearVariationalProblem> problem;

    }; // class DolfinPiece
    
} // namespace Modeling
} // namespace muq




#endif
