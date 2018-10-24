#include "MUQ/Optimization/NLoptOptimizer.h"

#include "MUQ/Utilities/Exceptions.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Optimization;


NLoptOptimizer::NLoptOptimizer(std::shared_ptr<CostFunction> cost,
                               pt::ptree const& pt) :
  OptimizerBase(cost, pt),
  algorithm(NLOptAlgorithm(pt.get<std::string>("Algorithm"))),
  opt(cost) {}


NLoptOptimizer::~NLoptOptimizer() {}


nlopt_algorithm NLoptOptimizer::NLOptAlgorithm(std::string const& alg) const {
  if( alg.compare("DIRECT")==0 ) { return NLOPT_GN_DIRECT; }
  if( alg.compare("DIRECTL")==0 ) { return NLOPT_GN_DIRECT_L; }
  if( alg.compare("CRS")==0 ) { return NLOPT_GN_CRS2_LM; }
  if( alg.compare("MLSL")==0 ) { return NLOPT_G_MLSL_LDS; }
  if( alg.compare("ISRES")==0 ) { return NLOPT_GN_ISRES; }
  if( alg.compare("COBYLA")==0 ) { return NLOPT_LN_COBYLA; }
  if( alg.compare("BOBYQA")==0 ) { return NLOPT_LN_BOBYQA; }
  if( alg.compare("NEWUOA")==0 ) { return NLOPT_LN_NEWUOA_BOUND; }
  if( alg.compare("PRAXIS")==0 ) { return NLOPT_LN_PRAXIS; }
  if( alg.compare("NM")==0 ) { return NLOPT_LN_NELDERMEAD; }
  if( alg.compare("SBPLX")==0 ) { return NLOPT_LN_SBPLX; }
  if( alg.compare("MMA")==0 ) { return NLOPT_LD_MMA; }
  if( alg.compare("SLSQP")==0 ) { return NLOPT_LD_SLSQP; }
  if( alg.compare("LBFGS")==0 ) { return NLOPT_LD_LBFGS; }
  if( alg.compare("PreTN")==0 ) { return NLOPT_LD_TNEWTON_PRECOND_RESTART; }
  if( alg.compare("LMVM")==0 ) { return NLOPT_LD_VAR2; }

  return NLOPT_LN_COBYLA;
}

void NLoptOptimizer::EvaluateImpl(ref_vector<boost::any> const& inputs) {

  // create the optimizer
  auto solver = nlopt_create(algorithm, opt->inputSizes(0));
  nlopt_set_min_objective(solver, Cost, &opt);

  if (ineqConstraints) {

    double ineqTol[ineqConstraints->numOutputs];
    Eigen::Map<const Eigen::VectorXd> ineqTolmap(ineqTol, ineqConstraints->numOutputs);
    const Eigen::VectorXd ineqTolEig = constraint_tol*Eigen::VectorXd::Ones(ineqConstraints->numOutputs);

    nlopt_add_inequality_mconstraint(solver,
                                     ineqConstraints->numOutputs,
                                     Constraint,
                                     &ineqConstraints,
                                     ineqTol); 
    
  }

  if (eqConstraints) {
    
    double eqTol[eqConstraints->numOutputs];
    Eigen::Map<const Eigen::VectorXd> eqTolmap(eqTol, eqConstraints->numOutputs);
    const Eigen::VectorXd eqTolEig = constraint_tol*Eigen::VectorXd::Ones(eqConstraints->numOutputs);

    nlopt_add_equality_mconstraint(solver,
                                   eqConstraints->numOutputs,
                                   Constraint,
                                   &eqConstraints,
                                   eqTol);

  }

  // set the tolerances
  nlopt_set_ftol_rel(solver, ftol_rel);
  nlopt_set_ftol_abs(solver, ftol_abs);
  nlopt_set_xtol_rel(solver, xtol_rel);
  nlopt_set_xtol_abs1(solver, xtol_abs);
  nlopt_set_maxeval(solver, maxEvals);

  const Eigen::VectorXd& xinit = boost::any_cast<Eigen::VectorXd const&>(inputs.at(0));

  outputs.resize(2);
  outputs[0] = xinit;
  Eigen::VectorXd& xopt = boost::any_cast<Eigen::VectorXd&>(outputs.at(0));

  double minf;

  const nlopt_result check = nlopt_optimize(solver, xopt.data(), &minf);

  if( check<0 )
    muq::ExternalLibraryError("NLOPT has failed with flag " + std::to_string(check));

  outputs[1] = minf;
  
}
      
std::pair<Eigen::VectorXd, double>
NLoptOptimizer::Solve(std::vector<Eigen::VectorXd> const& input) {

  std::vector<boost::any> ins;
  for (auto i : input)
    ins.push_back(i);
    
  Evaluate(ins);

  return std::pair<Eigen::VectorXd, double>(boost::any_cast<Eigen::VectorXd const&>(outputs[0]),
                                            boost::any_cast<double const>(outputs[1]));

}

void NLoptOptimizer::AddInequalityConstraint(std::shared_ptr<ModPiece> const& ineq) {
  ineqConstraints = ineq;
  UpdateInputs(ineq->numInputs-1);
}

void NLoptOptimizer::AddEqualityConstraint(std::shared_ptr<ModPiece> const& eq) {
  eqConstraints = eq;
  UpdateInputs(eq->numInputs-1);
}

void NLoptOptimizer::UpdateInputs(unsigned int const numNewIns) {

  for( unsigned int i=numInputs; i<numInputs+numNewIns; ++i )
    inputTypes[i] = typeid(Eigen::VectorXd).name(); 

  numInputs += numNewIns;

}

double NLoptOptimizer::Cost(unsigned int n,
                            const double* x,
                            double* grad,
                            void* f_data) {

  // The constraint
  std::shared_ptr<CostFunction> opt = *((std::shared_ptr<CostFunction>*) f_data);

  Eigen::Map<const Eigen::VectorXd> xmap(x, n);
  const Eigen::VectorXd& xeig = xmap;
  
  if (grad) {

    Eigen::Map<Eigen::VectorXd> gradmap(grad, n);
    const Eigen::VectorXd& gradeig =
      opt->Gradient(0, xeig, (Eigen::VectorXd)Eigen::VectorXd::Ones(1));
    gradmap = gradeig;

  }

  return opt->Cost(xeig);

}


void NLoptOptimizer::Constraint(unsigned int m,
                                double* result,
                                unsigned int n,
                                const double* x,
                                double* grad,
                                void* f_data) {

  // The constraint
  std::shared_ptr<ModPiece> opt =
    *((std::shared_ptr<ModPiece>*) f_data);

  Eigen::Map<const Eigen::VectorXd> xmap(x, n);
  const Eigen::VectorXd& xeig = xmap;

  if( grad ) {

    Eigen::Map<Eigen::MatrixXd> gradmap(grad, n, m);
    const Eigen::MatrixXd& gradeig =
      opt->Jacobian(0, 0, xeig);

    gradmap = gradeig;

  }

  Eigen::Map<Eigen::VectorXd> resultmap(result, m);
  std::vector<Eigen::VectorXd> resulteig = opt->Evaluate(xeig);
  resultmap = resulteig.at(0);

}

