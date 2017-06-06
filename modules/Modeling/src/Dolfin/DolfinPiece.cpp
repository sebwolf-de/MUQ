#include "MUQ/Modeling/Dolfin/DolfinPiece.h"


DolfinPiece::DolfinPiece(pybind11::object const& problemIn,
                         pybind11::object const& outputFieldIn,
                         std::vector<pybind11::object> inputs)
{
    ExtractInputs(inputs);
            
    // First, make sure the input object is actually an instance of LinearVariationalProblem
    CheckProblemType(problemIn, "LinearVariationalProblem");
    CheckProblemType(outputFieldIn, "Function");
            
    problem = SwigExtract(problemIn);
    outputField = SwigExtract(outputFieldIn);
};


void DolfinPiece::EvaluateFunc(std::vector<pybind11::object> const& inputs)
{
    // Set all of the expressions
    for(auto& expr : inputExprs)
    {
        CheckProblemType(inputs.at(expr.first),"list");
            
        auto xList = pybind11::list(inputs.at(expr.first));
        for(unsigned i=0; i<expr.second.second.size(); ++i)
            expr.second.first.attr(expr.second.second.at(i).c_str())  = xList[i];
    }

    // Set all the functions
    for(auto& f : inputFuncs)
    {
        CheckProblemType(inputs.at(f.first),"Function");
        *f.second = *SwigExtract(inputs.at(f.first)).Cast<std::shared_ptr<dolfin::Function>>();
    }
        
    dolfin::LinearVariationalSolver solver(problem);
    solver.solve();
}


Eigen::VectorXd DolfinPiece::EvaluateVec(Eigen::Ref<const Eigen::VectorXd> const& x, std::vector<double> vals)
{
    //Eigen::VectorXd& inVec = std::dynamic_pointer_cast<dolfin::EigenVector>(inputField->vector())->vec();
    //assert(inVec);

    // Update the value of the input field
    //inVec = x;
    //inputField->update();

    // Update each of the expressions
    //exprs.at(0).attr("bc_val") = leftVal;
    //exprs.at(1).attr("bc_val") = rightVal;
        
    dolfin::LinearVariationalSolver solver(problem);
    solver.solve();
        
    std::shared_ptr<dolfin::EigenVector> vec = std::dynamic_pointer_cast<dolfin::EigenVector>(outputField->vector());
    assert(vec);
        
    return vec->vec();
};


void DolfinPiece::CheckType(pybind11::object const& obj, std::string const& requiredType)
{
    std::string typeName = pybind11::handle(obj).ptr()->ob_type->tp_name;
        
    if(typeName.compare(requiredType))
    {
        throw std::invalid_argument("FenicsPiece constructor was given an instance of \"" + typeName + "\" but requires an instance of \"" + requiredType + "\"");
    }
}


void DolfinPiece::ExtractInputs(std::vector<pybind11::object> const& inputs){

    for(unsigned i=0; i<inputs.size(); ++ i){
            
        // If the input is a list or tuple, then it should define an expression
        if(pybind11::isinstance<pybind11::list>(inputs.at(i))){

            pybind11::list input(inputs[i]);
                 
            pybind11::object expr = input[0];
            CheckType(expr, "CompiledExpression");

            CheckType(input[1],"list");
            pybind11::list part2(input[1]);
                
            if(pybind11::isinstance<pybind11::list>(part2) || pybind11::isinstance<pybind11::tuple>(part2)){

                std::vector<std::string> names;
                for(auto& name : part2){
                    names.push_back(name.cast<std::string>());
                }
                inputExprs[i] = std::make_pair(expr, names);
            }
            
        }else{
            CheckType(inputs.at(i),"Function");
            inputFuncs[i] = SwigExtract(inputs.at(i));
        }
    }
}
