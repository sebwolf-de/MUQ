#ifndef SAMPLECOLLECTION_H
#define SAMPLECOLLECTION_H

#include "MUQ/SamplingAlgorithms/SamplingState.h"

#include <vector>
#include <Eigen/Core>

namespace muq{
  namespace SamplingAlgorithms{

    class SamplingStateIdentity{
    public:
      SamplingStateIdentity(int blockIndIn) : blockInd(blockIndIn){};

      Eigen::VectorXd const& operator()(SamplingState const& a);

      const int blockInd;

    private:
      Eigen::VectorXd output;
    };

    class SamplingStatePartialMoment{
    public:
      SamplingStatePartialMoment(int                    blockIndIn,
                                 int                    momentPowerIn,
                                 Eigen::VectorXd const& muIn) : blockInd(blockIndIn), momentPower(momentPowerIn), mu(muIn){};

      Eigen::VectorXd const& operator()(SamplingState const& a);

      const int blockInd;
      const int momentPower;
      const Eigen::VectorXd& mu;

    private:
      Eigen::VectorXd output;
    };

    /** @brief A class to hold and analyze a collection of SamplingState's
    */
    class SampleCollection{
    public:

      void Add(std::shared_ptr<SamplingState> newSamp);

      std::shared_ptr<SamplingState> at(unsigned i);
      const std::shared_ptr<SamplingState> at(unsigned i) const;

      //  Computes the componentwise central moments (e.g., variance, skewness, kurtosis, etc..) of a specific order
      boost::any CentralMoment(unsigned order, int blockDim=-1) const;

      boost::any Mean(int blockDim=-1) const;
      boost::any Variance(int blockDim=-1) const{return CentralMoment(2,blockDim);};
      boost::any Covariance(int blockDim=-1) const;

      /** @brief Returns the effective sample size of the samples
          @details For almost all random variables of interest, the central limit
                   theorem states that the variance of a Monte Carlo estimator
                   \f[\hat{\mu}_N = \frac{1}{N}\sum_{i=1}^N \theta^{(i)},\f]
                   satisfies
                   \f[ \sqrt{N}\mathbb{V}\left[\hat{\mu}_N-\mu\right] \rightarrow N(0,\sigma^2),\f]
                   as the number of samples \f$N\f$ increases.  Here, \f$mu\f$ is
                   the true mean of the random variable, \f$\sigma^2\f$ is the true
                   variance of the random variable.  This assumes that each \f$\theta^{(i)}\f$
                   is an independent sample.

                   When the samples are not independent,
                   the variance of the estimator will generally be larger.   The
                   "effective sample size" (ESS) describes how many independent
                   samples would have been needed to obtain the same estimator
                   variance.  In particular, let \f$\tilde{\mu}_N\f$ be a Monte
                   Carlo estimator based on \f$N\f$ correlated samples.  The ESS
                   is then given by the ratio of the esimator variances:
                   \f[ ESS = N \frac{\mathbb{V}[\hat{\mu}-\mu]}{\mathbb{V}[\tilde{\mu}-\mu]}. \f]
                   This function returns an approximation of this ESS using one
                   of two approaches:
                   - If the samples are correlated (i.e., the come from an MCMC algorithm)
      */
      double ESS(bool correlated=false) const;

    private:

      std::vector<std::shared_ptr<SamplingState>> samples;

      template<typename FuncType>
      static std::pair<double,Eigen::VectorXd> RecursiveSum(std::vector<const std::shared_ptr<SamplingState>>::iterator                         start,
                                                            std::vector<const std::shared_ptr<SamplingState>>::iterator                         end,
                                                            FuncType& f)
      {
        int numSamps = std::distance(start,end);
        const int maxSamps = 20;

        // If the number of samples is small enough, we can safely add them up directly
        if(numSamps<maxSamps){

          Eigen::VectorXd sum = (*start)->weight * f(**start);
          double weightSum = (*start)->weight;

          for(auto it=start+1; it!=end; ++it){
              sum += (*it)->weight * f(**it);
              weightSum += (*it)->weight;
          }
          return std::make_pair(weightSum, sum);

        // Otherwise, it's more numerically stable to add things pairwise
        }else{
          int halfDist = std::floor(0.5*numSamps);
          double weight1, weight2;
          Eigen::VectorXd sum1, sum2;
          std::tie(weight1,sum1) = RecursiveSum(start, start+halfDist, f);
          std::tie(weight2,sum2) = RecursiveSum(start+halfDist, end, f);

          return std::make_pair(weight1+weight2, (sum1+sum2).eval());
        }
      }


    };

  }
}



#endif // SAMPLECOLLECTION_H
