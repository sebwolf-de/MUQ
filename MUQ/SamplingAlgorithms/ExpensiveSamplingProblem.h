#ifndef EXPENSIVESAMPLINGPROBLEM_H_
#define EXPENSIVESAMPLINGPROBLEM_H_

#include "MUQ/Approximation/Regression/LocalRegression.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"

namespace muq {
  namespace SamplingAlgorithms {
    class ExpensiveSamplingProblem : public SamplingProblem {
    public:

      ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target, boost::property_tree::ptree& pt);

      ~ExpensiveSamplingProblem() = default;

      virtual double LogDensity(unsigned int const t, std::shared_ptr<SamplingState> state) override;

      unsigned int CacheSize() const;
      
    private:

      /**
	 @param[in] step The current MCMC step
	 @param[in] state The state where we are evalauting the log target
	 @param[out] neighbors The nearest neighbors
	 @param[out] results The log-target at the nearest neighbors
       */
      void RefineSurrogate(unsigned int const step, std::shared_ptr<SamplingState> state, std::vector<Eigen::VectorXd>& neighbors, std::vector<Eigen::VectorXd>& results);

      /**
	 Replace the point that is farthest from the new point
	 @param[in] state The point where we are evalauting the log target
	 @param[in] index The index of the point that we are replacing
	 @param[out] neighbors The nearest neighbors
	 @param[out] results The log-target at the nearest neighbors
       */
      void RefineSurrogate(Eigen::VectorXd const& point, unsigned int const index, std::vector<Eigen::VectorXd>& neighbors, std::vector<Eigen::VectorXd>& results) const;

      std::shared_ptr<muq::Approximation::LocalRegression> reg;

      /// Parameters for random refinement
      /**
	 Refine with probability \f$\beta = \beta_0 t^{-\beta_1}\f$.  \f$\beta_0\f$ is "BetaScale" and it defaults to \f$0\f$.  \f$\beta_1\f$ is "BetaExponent" and it defaults to \f$RAND_MAX\f$.
       */
      std::pair<double, double> beta;

      /// Level scaling for sturctural error
      double phi;

      /// The upper bound for the poisedness constant
      double lambda;

      /// Parameters for structural refinement
      /**
	 Refine if the error threshold exceeds \f$\gamma = \gamma_0 l^{-\gamma_1}\f$, where \f$l\f$ is the current error threshold level (ExpensiveSamplingProblem::level).  \f$\gamma_0\f$ is "GammaScale" and it defaults to \f$1\f$.  \f$\gamma_1\f$ is "GammaExponent" and it defaults to \f$1.0\f$.
       */
      std::pair<double, double> gamma;

      /// The current error threshold level
      unsigned int level = 1;
      
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
