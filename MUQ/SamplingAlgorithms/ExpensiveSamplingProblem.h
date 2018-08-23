#ifndef EXPENSIVESAMPLINGPROBLEM_H_
#define EXPENSIVESAMPLINGPROBLEM_H_

#include "MUQ/Approximation/Regression/LocalRegression.h"

#include "MUQ/SamplingAlgorithms/SamplingProblem.h"

namespace muq {
  namespace SamplingAlgorithms {
    class ExpensiveSamplingProblem : public SamplingProblem {
    public:

      ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target, boost::property_tree::ptree const& pt);

      ~ExpensiveSamplingProblem() = default;

      virtual double LogDensity(unsigned int const t, std::shared_ptr<SamplingState> state) override;
      
    private:

      void RefineSurrogate();

      std::shared_ptr<muq::Approximation::LocalRegression> reg;

      /// Parameters for random refinement
      /**
	 Refine with probability \f$\beta = \beta_0 t^{-\beta_1}\f$.  \f$\beta_0\f$ is "BetaScale" and it defaults to \f$0\f$.  \f$\beta_1\f$ is "BetaExponent" and it defaults to \f$RAND_MAX\f$.
       */
      std::pair<double, double> beta;

      /// Level scaling for sturctural error
      double phi1;

      /// Parameters for structural refinement
      /**
	 Refine if the error threshold exceeds \f$\gamma = \gamma_0 l^{\gamma_1}\f$, where \f$l\f$ is the current error threshold level (ExpensiveSamplingProblem::level).
       */
      std::pair<double, double> gamma;

      /// The current error threshold level
      unsigned int level = 1;
      
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
