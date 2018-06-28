#ifndef MarkovChain_H
#define MarkovChain_H

#include "MUQ/SamplingAlgorithms/SampleCollection.h"

namespace muq {
  namespace SamplingAlgorithms{

    /**
    @ingroup MCMC
    @class MarkovChain
    @brief A class for storing and working with the results of Markov chain Monte Carlo algorithms.
    @details The MarkovChain class is a child of SampleCollection where the sample
    weights correspond to the number of consecutive steps taking the same value,
    and the weights are unnormalized (i.e., do not sum to one).  This is a useful
    class for storing the chain produced by an MCMC algorithm without storing the
    duplicate points that result from rejected proposals.
    */
    class MarkovChain : public SampleCollection
    {
    public:

      virtual std::shared_ptr<SamplingState> at(unsigned i) override;
      virtual const std::shared_ptr<SamplingState> at(unsigned i) const override;

      virtual unsigned size() const override;

      /** Computes the effective sample size using the method described
      Ulli Wolff's "Monte Carlo errors with less error"
      */
      virtual Eigen::VectorXd ESS(int blockDim=-1) const override;

      /**
      Returns points in the Markov chain as a matrix.  Each column of the matrix
      is a different point.
      */
      virtual Eigen::MatrixXd AsMatrix(int blockDim=-1) const override;

      /** Returns a constant vector with valued \f$1/N\f$, where \f$N\f$ is the number
      of steps in the Markov chain.
      */
      virtual Eigen::VectorXd Weights() const override;

      static double SingleComponentESS(Eigen::Ref<const Eigen::VectorXd> const& trace);

    }; // class MarkovChain
  }
}

#endif
