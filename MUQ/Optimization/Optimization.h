#ifndef OPTIMIZATION_H_
#define OPTIMIZATION_H_

#include <boost/property_tree/ptree.hpp>

#include <nlopt.h>

#include "MUQ/Optimization/CostFunction.h"

namespace muq {
  namespace Optimization {
    class Optimization : public muq::Modeling::WorkPiece {
    public:

      Optimization(std::shared_ptr<CostFunction> cost, boost::property_tree::ptree const& pt);

      virtual ~Optimization();

      std::pair<Eigen::VectorXd, double> Solve(muq::Modeling::ref_vector<boost::any> const& inputs);

      template<typename ...Args>
	inline std::pair<Eigen::VectorXd, double> Solve(Args... args) {
	Evaluate(args...);

	return std::pair<Eigen::VectorXd, double>(boost::any_cast<Eigen::VectorXd const&>(outputs[0]), boost::any_cast<double const>(outputs[1]));
      }

    private:

      virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      static double Cost(unsigned int n, const double* x, double* grad, void* f_data);

      nlopt_algorithm NLOptAlgorithm(std::string const& alg) const;

      const nlopt_algorithm algorithm;

      struct OptHelper {
	OptHelper(std::shared_ptr<CostFunction> cost);

	virtual ~OptHelper();
	
	/// The cost function that we are trying to minimize
	std::shared_ptr<CostFunction> cost;

	muq::Modeling::ref_vector<Eigen::VectorXd> inputs;
      };

      OptHelper opt;

      const double ftol_rel, ftol_abs, xtol_rel, xtol_abs;
      const unsigned int maxEvals;
    };
  } // namespace Optimization
} // namespace muq

#endif
