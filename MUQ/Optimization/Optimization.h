#ifndef OPTIMIZATION_H_
#define OPTIMIZATION_H_

#include <boost/property_tree/ptree.hpp>

#include "MUQ/Optimization/CostFunction.h"

namespace muq {
namespace Optimization {
  /// Solve an optimization problem
  /**
     \f{eqnarray}{
     c &=& \min{J(x; \theta_1, ..., \theta_1)} \        \
     f_i(x) &\leq& 0 \                                  \
     g_i(x) &=& 0
     \f}
  */
  class Optimization : public muq::Modeling::WorkPiece {
  public:

    Optimization(std::shared_ptr<CostFunction> cost,
                 boost::property_tree::ptree const& pt);

    virtual ~Optimization();

    /// Add an inequality constraint to the optimization
    /**
       @param[in] ineq The constraint
    */
    virtual void AddInequalityConstraint(std::shared_ptr<CostFunction> ineq);
    
    /// Add an equality constraint to the optimization
    /**
       NOTE: the NLOPT algorithm used must be able to handle equality constraints
       @param[in] ineq The constraint
    */
    virtual void AddEqualityConstraint(std::shared_ptr<CostFunction> eq);

    
    /// Solve the optimization problem
    /**
       @param[in] inputs The first input is the variable we are optimizing over, then inputs to the cost function, and inputs to the constraints in the order they were added
       \return First: the argmin, second: the minimum cost
    */
    virtual std::pair<Eigen::VectorXd, double>
    Solve(muq::Modeling::ref_vector<boost::any> const& inputs)=0;

    /// Solve the optimization problem
    /**
       @param[in] args The first input is the variable we are optimizing over, then inputs to the cost function, and inputs to the constraints in the order they were added
       \return First: the argmin, second: the minimum cost
    */
    template<typename ...Args>
    inline std::pair<Eigen::VectorXd, double> Solve(Args... args) {

      Evaluate(args...);

      return std::pair<Eigen::VectorXd, double>(boost::any_cast<Eigen::VectorXd const&>(outputs[0]),
                                                boost::any_cast<double const>(outputs[1]));

    }

  protected:

    /// Update the inputs if a constraint is added
    /**
       Adding a constraint (potentially) increases the number of inputs to the optimization problem.  If the constraint requires inputs, add them to the optimization.
       @param[in] numNewIns The number of inputs (not the state) that the constraint requires
    */
    virtual void UpdateInputs(unsigned int const numNewIns)=0;
    
    
    /// Relative and absoluste tolerances on the cost function value and on the difference between successive values of the state
    const double ftol_rel, ftol_abs, xtol_rel, xtol_abs;
    
    /// Tolerance on the constraints
    const double constraint_tol;
    
    /// Maximum number of cost function evaluations
    const unsigned int maxEvals;
  };

} // namespace Optimization
} // namespace muq

#endif
