#ifndef ABSTRACTSAMPLINGPROBLEM_H
#define ABSTRACTSAMPLINGPROBLEM_H

#include <vector>
#include <assert.h>

#include <memory>

#include <Eigen/Core>

namespace muq{
  namespace SamplingAlgorithms{

    class SamplingState;

    /** @brief Abstract base class for MCMC and Importance Sampling problems.
        @details Assumes we are interested in sampling some probability density
        \f$\pi(x_1, x_2, \ldots, x_N)\f$, where each input \f$x_i\f$ is a vector
        valued block.  For instance, \f$x_1\f$ could be a Gaussian random variable
        and \f$x_2\f$ could be a vector containing the mean and variance hyperparameters.
        The block structure is primarily used in Gibbs or Metropolis-in-Gibbs
        MCMC algorithms.

        Both the number of blocks and the size of the blocks are stored.  However,
        if one of the blocks is not a vector type, it is safe to set the block size to -1.
     */
    class AbstractSamplingProblem
    {
    public:
      AbstractSamplingProblem(Eigen::VectorXi const& blockSizesIn) : numBlocks(blockSizesIn.size()),
                                                                     blockSizes(blockSizesIn)
                                                                     {assert(blockSizes.size()==numBlocks);};

      virtual ~AbstractSamplingProblem() = default;

      virtual double LogDensity(std::shared_ptr<SamplingState> state) = 0;

      virtual Eigen::VectorXd GradLogDensity(std::shared_ptr<SamplingState> state,
                                             unsigned                       blockWrt) = 0;

      const int numBlocks;
      const Eigen::VectorXi blockSizes;
    };

  }
}


#endif // ifndef ABSTRACTSAMPLINGPROBLEM_H
