#include "MUQ/SamplingAlgorithms/SampleCollection.h"
#include "MUQ/Utilities/AnyHelpers.h"

using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

Eigen::VectorXd const& SamplingStateIdentity::operator()(SamplingState const& a)
{
  if(blockInd<0){
    const int totalSize = a.TotalDim();
    const int numBlocks = a.state.size();

    output.resize(totalSize);
    int currInd = 0;
    for(int i=0; i<numBlocks; ++i){
      Eigen::VectorXd const& aVec = AnyConstCast(a.state.at(i));
      output.segment(currInd,aVec.size()) = aVec;
      currInd += aVec.size();
    }
    return output;

  }else{
    output.resize(0);
    Eigen::VectorXd const& aVec = AnyConstCast(a.state.at(blockInd));
    return aVec;
  }
}

Eigen::VectorXd const& SamplingStatePartialMoment::operator()(SamplingState const& a)
{
  if(blockInd<0){
    const int totalSize = a.TotalDim();
    const int numBlocks = a.state.size();

    output.resize(totalSize);
    int currInd = 0;
    for(int i=0; i<numBlocks; ++i){
      Eigen::VectorXd const& aVec = AnyConstCast(a.state.at(i));
      output.segment(currInd,aVec.size()) = (aVec-mu.segment(currInd,aVec.size())).array().pow(momentPower).matrix();
      currInd += aVec.size();
    }
    return output;

  }else{
    Eigen::VectorXd const& aVec = AnyConstCast(a.state.at(blockInd));
    output = (aVec-mu).array().pow(momentPower).matrix();
    return output;
  }
}

void SampleCollection::Add(std::shared_ptr<SamplingState> newSamp)
{
  samples.push_back(newSamp);
}

std::shared_ptr<SamplingState> SampleCollection::at(unsigned i)
{
  return samples.at(i);
}
const std::shared_ptr<SamplingState> SampleCollection::at(unsigned i) const
{
  return samples.at(i);
}

//  Computes the componentwise central moments (e.g., variance, skewness, kurtosis, etc..) of a specific order
boost::any SampleCollection::CentralMoment(unsigned order, int blockNum) const
{
  boost::any muAny = Mean(blockNum);
  Eigen::VectorXd const& mu = AnyConstCast(muAny);
  SamplingStatePartialMoment op(blockNum, order, mu);

  Eigen::VectorXd stateSum;
  double weightSum;

  std::tie(weightSum, stateSum) = RecursiveSum(samples.begin(), samples.end(), op);
  return (stateSum / weightSum).eval();
}

boost::any SampleCollection::Mean(int blockNum) const
{
    SamplingStateIdentity op(blockNum);

    Eigen::VectorXd stateSum;
    double weightSum;


    std::tie(weightSum, stateSum) = RecursiveSum(samples.begin(), samples.end(), op);
    return (stateSum / weightSum).eval();
}

boost::any SampleCollection::Covariance(int blockInd) const
{
  const int numSamps = samples.size();
  Eigen::MatrixXd samps;
  Eigen::VectorXd weights(numSamps);

  if(blockInd<0){

    const int totalSize = samples.at(0)->TotalDim();
    const int numBlocks = samples.at(0)->state.size();

    samps.resize(totalSize, numSamps);

    for(int i=0; i<numSamps; ++i){
      weights(i) = samples.at(i)->weight;

      int currInd = 0;
      for(int block = 0; block<numBlocks; ++block){
        Eigen::VectorXd const& tempVec  = AnyConstCast(samples.at(i)->state.at(block));
        samps.col(i).segment(currInd, tempVec.size()) = tempVec;
        currInd += tempVec.size();
      }

    }

  }else{


    const int blockSize = boost::any_cast<Eigen::VectorXd const&>(samples.at(0)->state.at(blockInd)).size();

    samps.resize(blockSize, numSamps);

    for(int i=0; i<numSamps; ++i){
      weights(i) = samples.at(i)->weight;
      samps.col(i) = boost::any_cast<Eigen::VectorXd const&>(samples.at(i)->state.at(blockInd));
    }
  }

  Eigen::VectorXd mu = (samps * weights) / weights.sum();
  Eigen::MatrixXd cov = (samps.colwise() - mu) * weights.asDiagonal() * (samps.colwise()-mu).transpose() / weights.sum();
  return cov;
}

// std::pair<double,Eigen::VectorXd> SampleCollection::RecursiveSum(std::vector<const std::shared_ptr<SamplingState>>::iterator       start,
//                                                                  std::vector<const std::shared_ptr<SamplingState>>::iterator       end,
//                                                                  std::function<Eigen::VectorXd(boost::any)> f)
// {
//   int numSamps = std::distance(start,end);
//   const int maxSamps = 20;
//
//   // If the number of samples is small enough, we can safely add them up directly
//   if(numSamps<maxSamps){
//
//     Eigen::VectorXd sum = (*start)->weight * f(**start);
//     double weightSum = (*start)->weight;
//
//     for(auto it=start+1; it!=end; ++it){
//         sum += (*it)->weight * f(**it);
//         weightSum += (*it)->weight;
//     }
//     return std::make_pair(weightSum, sum);
//
//   // Otherwise, it's more numerically stable to add things pairwise
//   }else{
//     int halfDist = std::floor(0.5*numSamps);
//     double weight1, weight2;
//     Eigen::VectorXd sum1, sum2;
//     std::tie(weight1,sum1) = RecursiveSum(start, start+halfDist);
//     std::tie(weight2,sum2) = RecursiveSum(start+halfDist, end);
//
//     return std::make_pair(weight1+weight2, (sum1+sum2).eval());
//   }
// }
