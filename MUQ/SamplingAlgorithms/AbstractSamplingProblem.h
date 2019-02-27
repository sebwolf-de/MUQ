#ifndef ABSTRACTSAMPLINGPROBLEM_H
#define ABSTRACTSAMPLINGPROBLEM_H

#include <vector>
#include <assert.h>

#include <memory>

#include <Eigen/Core>
#include <boost/property_tree/ptree.hpp>

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

      enum SampleType {
	       Proposed,
	        Accepted
      };

      inline AbstractSamplingProblem(Eigen::VectorXi const& blockSizesIn) : numBlocks(blockSizesIn.size()),
                                                                     blockSizes(blockSizesIn)
                                                                     {assert(blockSizes.size()==numBlocks);};

      virtual ~AbstractSamplingProblem() = default;

      virtual double LogDensity(unsigned int const t, std::shared_ptr<SamplingState> state, AbstractSamplingProblem::SampleType type) = 0;

      /** Sometimes, there will be problem-specific options that need to be passed
          to the SamplingAlgorithm.  This function adds any of those options to the
          given property_tree.
      */
      virtual void AddOptions(boost::property_tree::ptree & pt) const{};

      const int numBlocks;
      const Eigen::VectorXi blockSizes;

    };

  }
}


#endif // ifndef ABSTRACTSAMPLINGPROBLEM_H
