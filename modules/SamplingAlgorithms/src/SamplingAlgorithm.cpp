#include "MUQ/SamplingAlgorithms/SamplingAlgorithm.h"

#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

SamplingAlgorithm::SamplingAlgorithm() : WorkPiece(std::map<unsigned int, std::string>({
      std::pair<unsigned int, std::string>(0, typeid(pt::ptree).name()), // algorithm parameters
	std::pair<unsigned int, std::string>(1, typeid(std::shared_ptr<SamplingProblem>).name()) // sampling problem
	})) {}

SamplingAlgorithm::~SamplingAlgorithm() {}

void SamplingAlgorithm::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // get the options for this algorithm
  pt::ptree pt = boost::any_cast<pt::ptree>(inputs[0]);

  // get the problem that defines the distribution we are trying to characterize/sample
  std::shared_ptr<SamplingProblem> problem = boost::any_cast<std::shared_ptr<SamplingProblem> >(inputs[1]);
  assert(problem);

  // get the kernel that generates the next sample
  std::shared_ptr<TransitionKernel> kernel = Kernel(pt, problem);

  // get the number of samples
  const unsigned int T = pt.get<unsigned int>("SamplingAlgorithm.NumSamples");

  for( unsigned int t=0; t<T; ++t ) { // loop through each sample
    //SampleOnce(problem);
  }
}

/*boost::any SamplingAlgorithm::SampleOnce(std::shared_ptr<SamplingProblem> problem) const {
  problem->Sample();
  
  return boost::none;
  }*/

