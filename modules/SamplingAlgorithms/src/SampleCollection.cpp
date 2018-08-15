#include "MUQ/SamplingAlgorithms/SampleCollection.h"

#include "MUQ/Utilities/HDF5/HDF5File.h"
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

std::pair<double,double> SampleCollection::RecursiveWeightSum(std::vector<std::shared_ptr<SamplingState>>::const_iterator start,
                                                              std::vector<std::shared_ptr<SamplingState>>::const_iterator end)
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

void SampleCollection::WriteToFile(std::string const& filename, std::string const& dataset) const {
  // open the hdf5 file
  auto hdf5file = std::make_shared<HDF5File>(filename);

  // write the sample matrix and weights
  hdf5file->WriteMatrix(dataset+"/samples", AsMatrix());
  hdf5file->WriteMatrix(dataset+"/weights", (Eigen::RowVectorXd)Weights());

  // meta data
  std::unordered_map<std::string, Eigen::MatrixXd> meta;

  //for( auto samp : samples ) {
  for( unsigned int i=0; i<samples.size(); ++i ) {
    for( auto it = samples[i]->meta.begin(); it!=samples[i]->meta.end(); ++it ) {
      if( it->second.type()==typeid(Eigen::Vector2d) ) {
	// create a matrix for the meta data 
	if( meta.find(it->first)==meta.end() ) { meta.insert({it->first, Eigen::MatrixXd::Constant(2, samples.size(), std::numeric_limits<double>::quiet_NaN())}); }
	meta[it->first].col(i) = boost::any_cast<Eigen::Vector2d const&>(it->second);

	continue;		
      }
      
      if( it->second.type()==typeid(Eigen::Vector3d) ) {
	// create a matrix for the meta data 
	if( meta.find(it->first)==meta.end() ) { meta.insert({it->first, Eigen::MatrixXd::Constant(3, samples.size(), std::numeric_limits<double>::quiet_NaN())}); }
	meta[it->first].col(i) = boost::any_cast<Eigen::Vector3d const&>(it->second);

	continue;		
      }

      if( it->second.type()==typeid(Eigen::Vector4d) ) {
	// create a matrix for the meta data 
	if( meta.find(it->first)==meta.end() ) { meta.insert({it->first, Eigen::MatrixXd::Constant(4, samples.size(), std::numeric_limits<double>::quiet_NaN())}); }
	meta[it->first].col(i) = boost::any_cast<Eigen::Vector4d const&>(it->second);

	continue;		
      }

      if( it->second.type()==typeid(Eigen::VectorXd) ) {
	// create a matrix for the meta data 
	if( meta.find(it->first)==meta.end() ) { meta.insert({it->first, Eigen::MatrixXd::Constant(boost::any_cast<Eigen::VectorXd const&>(it->second).size(), samples.size(), std::numeric_limits<double>::quiet_NaN())}); }
	
	meta[it->first].col(i) = boost::any_cast<Eigen::VectorXd const&>(it->second);

	continue;		
      }
      
      // create a matrix, assuming scalar type, for the meta data
      if( meta.find(it->first)==meta.end() ) { meta.insert({it->first, Eigen::MatrixXd::Constant(1, samples.size(), std::numeric_limits<double>::quiet_NaN())}); }
      
      if( it->second.type()==typeid(double) ) { // doubles
	meta[it->first](i) = boost::any_cast<double const>(it->second);
	continue;
      }

      if( it->second.type()==typeid(float) ) { // floats
	meta[it->first](i) = boost::any_cast<float const>(it->second);
	continue;
      }

      if( it->second.type()==typeid(int) ) { // ints
	meta[it->first](i) = boost::any_cast<int const>(it->second);
	continue;
      }

      if( it->second.type()==typeid(unsigned int) ) { // unsigned ints
	meta[it->first](i) = boost::any_cast<unsigned int const>(it->second);
	continue;
      }
    }
  }

  // write meta data to file
  for( const auto& data : meta ) { hdf5file->WriteMatrix(dataset+"/"+data.first, data.second); }
}
