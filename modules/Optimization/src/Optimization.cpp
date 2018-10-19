#include "MUQ/Optimization/Optimization.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Optimization;


Optimization::Optimization(std::shared_ptr<CostFunction> cost,
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


Optimization::~Optimization() {}
