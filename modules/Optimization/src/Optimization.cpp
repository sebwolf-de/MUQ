#include "MUQ/Optimization/Optimization.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Optimization;


Optimization::Optimization(std::shared_ptr<CostFunction> cost,
                           pt::ptree const& pt) :
  WorkPiece(cost->InputTypes(),
            cost->numInputs,
            std::vector<std::string>({typeid(Eigen::VectorXd).name(), typeid(double).name()})),
  opt(cost, 1),  
  ftol_rel(pt.get<double>("Ftol.AbsoluteTolerance", 1.0e-8)),
  ftol_abs(pt.get<double>("Ftol.RelativeTolerance", 1.0e-8)),
  xtol_rel(pt.get<double>("Xtol.AbsoluteTolerance", 0.0)),
  xtol_abs(pt.get<double>("Xtol.RelativeTolerance", 0.0)),
  constraint_tol(pt.get<double>("ConstraintTolerance", 1.0e-8)),
  maxEvals(pt.get<unsigned int>("MaxEvaluations", 100)) {}


Optimization::~Optimization() {}

double Optimization::Cost(unsigned int n,
                          const double* x,
                          double* grad,
                          void* f_data) {

  
  CostHelper* opt = (CostHelper*)f_data;

  Eigen::Map<const Eigen::VectorXd> xmap(x, n);
  const Eigen::VectorXd& xeig = xmap;
  opt->inputs.at(0) = std::cref(xeig);
  
  if( grad ) {
    Eigen::Map<Eigen::VectorXd> gradmap(grad, n);
    const Eigen::VectorXd& gradeig =
      opt->cost->Gradient(0,
                          opt->inputs,
                          (Eigen::VectorXd)Eigen::VectorXd::Ones(1));
    gradmap = gradeig;
  }
  
  return opt->cost->Cost(opt->inputs);

}


Optimization::CostHelper::CostHelper(std::shared_ptr<CostFunction> cost,
                                     unsigned int const firstin) :
  cost(cost), firstin(firstin) {}

Optimization::CostHelper::~CostHelper() {}

void Optimization::CostHelper::SetInputs(ref_vector<boost::any> const& ins) {
  inputs.clear();
  
  // set initial condition
  inputs.push_back(std::cref(boost::any_cast<Eigen::VectorXd const&>(ins[0])));

  // set other inputs
  for( unsigned int i=firstin; i<firstin+cost->numInputs-1; ++i ) {
    inputs.push_back(std::cref(boost::any_cast<Eigen::VectorXd const&>(ins[i])));
  }
}
