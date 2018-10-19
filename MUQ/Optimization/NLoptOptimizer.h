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

      /// Add an inequality constraint to the optimization
      /**
         @param[in] ineq The constraint
      */
      void AddInequalityConstraint(std::shared_ptr<CostFunction> ineq);
      
      /// Add an equality constraint to the optimization
      /**
         NOTE: the NLOPT algorithm used must be able to handle equality constraints
         @param[in] ineq The constraint
      */
      void AddEqualityConstraint(std::shared_ptr<CostFunction> eq);

      virtual std::pair<Eigen::VectorXd, double>
      Solve(muq::Modeling::ref_vector<boost::any> const& inputs) override;

    private:

      
      /// Evaluate either the cost function or a constraint
      /**
         @param[in] n The size of the input
         @param[in] x The current point
         @param[out] grad The gradient of the cost/constraint 
         @param[in] f_data An Optimization::CostHelper
         \return The cost/constraint value
      */
      static double Cost(unsigned int n,
                         const double* x,
                         double* grad,
                         void* f_data);
      
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


      /// A structure to help evaluate the cost function and constraints
      struct CostHelper {
        /**
           @param[in] cost The muq::Optimization::CostFunction that evaluates either the cost function or the constraint
           @param[in] firstin The index of the optimziation inputs where this functions inputs begin
        */
        CostHelper(std::shared_ptr<CostFunction> cost, unsigned int const firstin);
        
        virtual ~CostHelper();
        
        /// Given an input to the optimization problem, set the inputs of this fucntion
        /**
           @param[in] ins The inputs to the optimization problem
        */
        void SetInputs(muq::Modeling::ref_vector<boost::any> const& ins);
        
        /// The cost function that we are trying to minimize or a cosntraint
        std::shared_ptr<CostFunction> cost;
        
        /// The index of the optimziation inputs where this functions inputs begin
        const unsigned int firstin;
        
	/// The inputs to this function
        muq::Modeling::ref_vector<Eigen::VectorXd> inputs;

      };

      /// The cost function that we are trying to minimize
      CostHelper opt;
      
      /// Inequality constraints
      std::vector<CostHelper> ineqConstraints;
      
      /// Equality constraints
      /**
         NOTE: the solver muq::Optimization::Optimization::algorithm must be able to handle equality constraints
      */
      std::vector<CostHelper> eqConstraints;
      
    }; // class NLoptOptimizer
      
  } // namespace Optimization
} // namespace muq

#endif
