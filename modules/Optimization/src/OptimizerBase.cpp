#include "MUQ/Optimization/OptimizerBase.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Optimization;


OptimizerBase::OptimizerBase(std::shared_ptr<CostFunction> cost,
                             pt::ptree const& pt) :
  WorkPiece(cost->InputTypes(),
            cost->numInputs,
            std::vector<std::string>({typeid(Eigen::VectorXd).name(), typeid(double).name()})),
  ftol_rel(pt.get<double>("Ftol.AbsoluteTolerance", 1.0e-8)),
  ftol_abs(pt.get<double>("Ftol.RelativeTolerance", 1.0e-8)),
  xtol_rel(pt.get<double>("Xtol.AbsoluteTolerance", 0.0)),
  xtol_abs(pt.get<double>("Xtol.RelativeTolerance", 0.0)),
  constraint_tol(pt.get<double>("ConstraintTolerance", 1.0e-8)),
  maxEvals(pt.get<unsigned int>("MaxEvaluations", 100)) {}


OptimizerBase::~OptimizerBase() {}


void OptimizerBase::AddInequalityConstraint(std::shared_ptr<ModPiece> const& ineq) {
  std::cerr << "ERROR: This optimization method does not handle inequality constraints." << std::endl;
}


void OptimizerBase::AddEqualityConstraint(std::shared_ptr<ModPiece> const& eq) {
  std::cerr << "ERROR: This optimization method does not handle equality constraints." << std::endl;
}
  
