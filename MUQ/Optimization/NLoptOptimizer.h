#ifndef NLOPTOPTIMIZATION_H_
#define NLOPTOPTIMIZATION_H_

#include <boost/property_tree/ptree.hpp>

#include <nlopt.h>

#include "MUQ/Optimization/Optimization.h"
#include "MUQ/Optimization/CostFunction.h"

namespace muq {
  namespace Optimization {

    class NLoptOptimizer : public Optimization {
    public:

      NLoptOptimizer(std::shared_ptr<CostFunction> cost,
                     boost::property_tree::ptree const& pt);

      virtual ~NLoptOptimizer();

      virtual void AddInequalityConstraint(std::shared_ptr<CostFunction> ineq) override;
      virtual void AddEqualityConstraint(std::shared_ptr<CostFunction> eq) override;

      virtual std::pair<Eigen::VectorXd, double>
      Solve(muq::Modeling::ref_vector<boost::any> const& inputs) override;

    private:

      /// Override the evaluate impl method (solve the optimization problem)
      /**
	 @param[in] args The first input is the variable we are optimizing over, then inputs to the cost function, and inputs to the constraints in the order they were added
       */
      virtual void
      EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      virtual void UpdateInputs(unsigned int const numNewIns) override;

      /// Get the NLOPT algorithm we are using
      /**
	 @param[in] alg User-input algorithm
	 \return The NLOPT algorithm
       */
      nlopt_algorithm NLOptAlgorithm(std::string const& alg) const;

      /// The algorithm used to solve the problem
      const nlopt_algorithm algorithm;

    }; // class NLoptOptimizer
      
  } // namespace Optimization
} // namespace muq

#endif
