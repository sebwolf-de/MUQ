#include "MUQ/SamplingAlgorithms/SampleCollection.h"

using namespace muq::SamplingAlgorithms;

void SampleCollection::Add(SamplingState const& newSamp)
{
  samples.push_back(newSamp);
}

SamplingState&       SampleCollection::at(unsigned i)
{
  return samples.at(i);
}

SamplingState const& SampleCollection::at(unsigned i) const
{
  return samples.at(i);
}

//  Computes the componentwise central moments (e.g., variance, skewness, kurtosis, etc..) of a specific order
boost::any SampleCollection::CentralMoment(unsigned order) const
{

}

boost::any SampleCollection::Mean() const
{

}

boost::any SampleCollection::Covariance() const
{

}
