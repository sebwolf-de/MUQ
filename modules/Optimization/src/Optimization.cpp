#include "MUQ/Optimization/Optimization.h"

using namespace muq::Modeling;
using namespace muq::Optimization;

Optimization::Optimization(std::shared_ptr<CostFunction> cost) :
  WorkPiece(cost->InputTypes(), cost->numInputs, std::vector<std::string>({typeid(Eigen::VectorXd).name(), typeid(double).name()})),
  opt(cost) {}

Optimization::~Optimization() {}

double Optimization::Cost(unsigned int n, const double* x, double* grad, void* f_data) {
  OptHelper* opt = (OptHelper*)f_data;

  Eigen::Map<const Eigen::VectorXd> xmap(x, n);
  const Eigen::VectorXd& xeig = xmap;
  //std::cout << xeig.transpose() << std::endl;
  opt->inputs.at(0) = std::cref(xeig);
  
  if( grad ) {
    Eigen::Map<Eigen::VectorXd> gradmap(grad, n);
    const Eigen::VectorXd& gradeig = opt->cost->Gradient(0, opt->inputs, (Eigen::VectorXd)Eigen::VectorXd::Ones(n));
    gradmap = gradeig;
  }
  
  return opt->cost->Cost(opt->inputs);
}

void Optimization::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  opt.inputs.clear();
  const Eigen::VectorXd& xinit = boost::any_cast<Eigen::VectorXd const&>(inputs.at(0));
  opt.inputs.push_back(std::cref(xinit));
  for( auto it=inputs.begin()+1; it!=inputs.end(); ++it ) {
    opt.inputs.push_back(std::cref(boost::any_cast<Eigen::VectorXd const&>(*it)));
  }
  
  // create the optimizer  
  auto solver = nlopt_create(NLOPT_LD_MMA, opt.cost->inputSizes(0));

  nlopt_set_min_objective(solver, Optimization::Cost, &opt);

  nlopt_set_ftol_rel(solver, 1e-14);
  nlopt_set_ftol_abs(solver, 1e-14);

  outputs.resize(2);
  outputs[0] = xinit;
  Eigen::VectorXd& xopt = boost::any_cast<Eigen::VectorXd&>(outputs.at(0));
  
  double minf;
  const double check = nlopt_optimize(solver, xopt.data(), &minf);
  assert(check>=0);
  
  outputs[1] = minf;
}
      
std::pair<Eigen::VectorXd, double> Optimization::Solve(ref_vector<boost::any> const& inputs) {
  Evaluate(inputs);

  return std::pair<Eigen::VectorXd, double>(boost::any_cast<Eigen::VectorXd const&>(outputs[0]), boost::any_cast<double const>(outputs[1]));
}

Optimization::OptHelper::OptHelper(std::shared_ptr<CostFunction> cost) : cost(cost) {}

Optimization::OptHelper::~OptHelper() {}
