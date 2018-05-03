#include "MUQ/SamplingAlgorithms/SamplingState.h"
#include <Eigen/Core>
#include "MUQ/Utilities/AnyHelpers.h"

using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

SamplingState::SamplingState(boost::any const& stateIn, double weight) : state({stateIn}), weight(weight) {}
SamplingState::SamplingState(std::vector<boost::any> const& stateIn, double const weight) : state(stateIn), weight(weight) {}

SamplingState::~SamplingState() {}

bool SamplingState::HasMeta(std::string const& metaKey)
{
  auto iter = meta.find(metaKey);
  return iter != meta.end();
}

int SamplingState::TotalDim() const
{
  int sum = 0;
  for(auto& s : state){
    Eigen::VectorXd const& temp = AnyConstCast(s);
    sum += temp.size();
  }
  return sum;
}
