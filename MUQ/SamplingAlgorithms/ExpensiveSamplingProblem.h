#ifndef EXPENSIVESAMPLINGPROBLEM_H_
#define EXPENSIVESAMPLINGPROBLEM_H_

#include "MUQ/config.h"

#if MUQ_HAS_PARCER
#include <parcer/Communicator.h>
#endif

#include "MUQ/Approximation/Regression/LocalRegression.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"

namespace muq {
  namespace SamplingAlgorithms {
    class ExpensiveSamplingProblem : public SamplingProblem {
    public:

      /**
	 @param[in] target The target distribution
	 @param[in] pt Options for the sampling problem
       */
      ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target, boost::property_tree::ptree& pt);

#if MUQ_HAS_PARCER
      /**
	 @param[in] target The target distribution
	 @param[in] pt Options for the sampling problem
	 @param[in] comm Parallel communicator for regression
       */
      ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target, boost::property_tree::ptree& pt, std::shared_ptr<parcer::Communicator> comm);
#endif

      ~ExpensiveSamplingProblem() = default;

      virtual double LogDensity(unsigned int const t, std::shared_ptr<SamplingState> state, AbstractSamplingProblem::SampleType type) override;

      unsigned int CacheSize() const;
      
    private:

      /// Set up the sampling problem
      /**
	 @param[in] pt Options for the sampling problem
       */
      void SetUp(boost::property_tree::ptree& pt);

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
      void RefineSurrogate(Eigen::VectorXd const& point, unsigned int const index, std::vector<Eigen::VectorXd>& neighbors, std::vector<Eigen::VectorXd>& results);

      /**
	 @param[in] state The point where we are evalauting the log target
       */
      void UpdateGlobalData(Eigen::VectorXd const& point);

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

      /// Exponenent for delta refinement
      double delta;

      /// Parameters for structural refinement
      /**
	 Refine if the error threshold exceeds \f$\gamma = \gamma_0 l^{-\gamma_1}\f$, where \f$l\f$ is the current error threshold level (ExpensiveSamplingProblem::level).  \f$\gamma_0\f$ is "GammaScale" and it defaults to \f$1\f$.  \f$\gamma_1\f$ is "GammaExponent" and it defaults to \f$1.0\f$.
       */
      std::pair<double, double> gamma;

      /// The current error threshold level
      unsigned int level = 1;

      /// Cumulative beta refinements
      unsigned int cumbeta = 0;

      /// Cumulative gamma refinements
      unsigned int cumgamma = 0;

      /// Cumulative kappa refinements
      unsigned int cumkappa = 0;
      
      /// Cumulative delta refinements
      unsigned int cumdelta = 0;

      /// Global mean of evaluated locations
      Eigen::VectorXd globalMean;

      double radius_avg = 0.0;

      /// Global radius of evaluated locations
      double radius_max = 0.0;
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
