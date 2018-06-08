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
      output.segment(currInd,a.state.at(i).size()) = a.state.at(i);
      currInd += a.state.at(i).size();
    }
    return output;

  }else{
    output.resize(0);
    return a.state.at(blockInd);
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
      output.segment(currInd,a.state.at(i).size()) = (a.state.at(i)-mu.segment(currInd,a.state.at(i).size())).array().pow(momentPower).matrix();
      currInd += a.state.at(i).size();
    }
    return output;

  }else{
    output = (a.state.at(blockInd)-mu).array().pow(momentPower).matrix();
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
Eigen::VectorXd SampleCollection::CentralMoment(unsigned order, int blockNum) const
{
  Eigen::VectorXd mu = Mean(blockNum);
  SamplingStatePartialMoment op(blockNum, order, mu);

  Eigen::VectorXd stateSum;
  double weightSum;

  std::tie(weightSum, stateSum) = RecursiveSum(samples.begin(), samples.end(), op);
  return (stateSum / weightSum).eval();
}

Eigen::VectorXd SampleCollection::Mean(int blockNum) const
{
    SamplingStateIdentity op(blockNum);

    Eigen::VectorXd stateSum;
    double weightSum;

    std::tie(weightSum, stateSum) = RecursiveSum(samples.begin(), samples.end(), op);
    return (stateSum / weightSum).eval();
}

Eigen::MatrixXd SampleCollection::Covariance(int blockInd) const
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
        samps.col(i).segment(currInd, samples.at(i)->state.at(block).size()) = samples.at(i)->state.at(block);
        currInd += samples.at(i)->state.at(block).size();
      }

    }

  }else{

    const int blockSize = samples.at(0)->state.at(blockInd).size();

    samps.resize(blockSize, numSamps);

    for(int i=0; i<numSamps; ++i){
      weights(i) = samples.at(i)->weight;
      samps.col(i) = samples.at(i)->state.at(blockInd);
    }
  }

  Eigen::VectorXd mu = (samps * weights) / weights.sum();
  Eigen::MatrixXd cov = (samps.colwise() - mu) * weights.asDiagonal() * (samps.colwise()-mu).transpose() / weights.sum();
  return cov;
}

std::pair<double,double> SampleCollection::RecursiveWeightSum(std::vector<const std::shared_ptr<SamplingState>>::iterator start,
                                                              std::vector<const std::shared_ptr<SamplingState>>::iterator end)
{
  int numSamps = std::distance(start,end);
  const int maxSamps = 20;

  // If the number of samples is small enough, we can safely add them up directly
  if(numSamps<maxSamps){

    double weightSquareSum = (*start)->weight * (*start)->weight;
    double weightSum = (*start)->weight;

    for(auto it=start+1; it!=end; ++it){
        weightSquareSum += (*it)->weight * (*it)->weight;
        weightSum += (*it)->weight;
    }
    return std::make_pair(weightSum, weightSquareSum);

  }else{
    int halfDist = std::floor(0.5*numSamps);
    double weight1, weight2, squaredSum1, squaredSum2;
    std::tie(weight1,squaredSum1) = RecursiveWeightSum(start, start+halfDist);
    std::tie(weight2,squaredSum2) = RecursiveWeightSum(start+halfDist, end);

    return std::make_pair(weight1+weight2, squaredSum1+squaredSum2);
  }
}


Eigen::VectorXd SampleCollection::ESS(int blockDim) const
{
  if(samples.size()==0)
    return Eigen::VectorXd();

  double weightSum = 0.0;
  double squaredSum = 0.0;
  std::tie(weightSum, squaredSum) = RecursiveWeightSum(samples.begin(), samples.end());

  int blockSize;
  if(blockDim<0){
    blockSize = samples.at(0)->TotalDim();
  }else{
    blockSize = samples.at(0)->state.at(blockDim).size();
  }
  return (weightSum*weightSum / squaredSum) * Eigen::VectorXd::Ones(blockSize);
}

Eigen::MatrixXd SampleCollection::AsMatrix(int blockDim) const
{
  if(samples.size()==0)
    return Eigen::MatrixXd();

  if(blockDim<0){

    Eigen::MatrixXd output(samples.at(0)->TotalDim(), samples.size());
    for(int i=0; i<samples.size(); ++i){
      int startInd = 0;

      for(int block=0; block<samples.at(0)->state.size(); ++block){
        int blockSize = samples.at(0)->state.at(block).size();
        output.col(i).segment(startInd, blockSize) = samples.at(i)->state.at(block);
        startInd += blockSize;
      }
    }

    return output;

  }else{
    Eigen::MatrixXd output(samples.at(0)->state.at(blockDim).size(), samples.size());
    for(int i=0; i<samples.size(); ++i)
      output.col(i) = samples.at(0)->state.at(blockDim);
    return output;
  }
}

Eigen::VectorXd SampleCollection::Weights() const
{
  Eigen::VectorXd output(samples.size());
  for(int i=0; i<samples.size(); ++i)
    output(i) = samples.at(i)->weight;
  return output;
}
