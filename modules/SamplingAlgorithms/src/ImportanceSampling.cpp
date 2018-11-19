#include "MUQ/SamplingAlgorithms/ImportanceSampling.h"

#include "MUQ/config.h"

#include "MUQ/Modeling/Distributions/Density.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;

ImportanceSampling::ImportanceSampling(std::shared_ptr<ModPiece> const& target, std::shared_ptr<Distribution> const& bias, pt::ptree const& pt) : SamplingAlgorithm(std::make_shared<SampleCollection>()),
  numSamps(pt.get<unsigned int>("NumSamples")), target(target), bias(bias) {}

ImportanceSampling::ImportanceSampling(std::shared_ptr<ModPiece> const& target, std::shared_ptr<Distribution> const& bias, std::vector<Eigen::VectorXd> hyperparameters, pt::ptree const& pt) : SamplingAlgorithm(std::make_shared<SampleCollection>()), numSamps(pt.get<unsigned int>("NumSamples")), target(target), bias(bias), hyperparameters(hyperparameters) {}

#if MUQ_HAS_PARCER
ImportanceSampling::ImportanceSampling(std::shared_ptr<ModPiece> const& target, std::shared_ptr<Distribution> const& bias, pt::ptree const& pt, std::shared_ptr<parcer::Communicator> const& comm) : SamplingAlgorithm(std::make_shared<SampleCollection>()),
  numSamps(pt.get<unsigned int>("NumSamples")), target(target), bias(bias), comm(comm) {}

ImportanceSampling::ImportanceSampling(std::shared_ptr<ModPiece> const& target, std::shared_ptr<Distribution> const& bias, std::vector<Eigen::VectorXd> hyperparameters, pt::ptree const& pt, std::shared_ptr<parcer::Communicator> const& comm) : SamplingAlgorithm(std::make_shared<SampleCollection>()), numSamps(pt.get<unsigned int>("NumSamples")), target(target), bias(bias), hyperparameters(hyperparameters), comm(comm) {}
#endif

ImportanceSampling::~ImportanceSampling() {}

void ImportanceSampling::RunSerial(std::vector<Eigen::VectorXd>& biasingPara) {
  // loop through the samples to generate them
  for( unsigned int i=0; i<numSamps; ++i ) {
    // sample from the biasing distribution
    biasingPara[0] = bias->Sample(hyperparameters);

    // compute the weight
    const double logweight = target->Evaluate(biasingPara[0])[0](0) - bias->LogDensity(biasingPara);

    // store the sample
    samples->Add(std::make_shared<SamplingState>(biasingPara[0], std::exp(logweight)));
  }
}

#if MUQ_HAS_PARCER
struct ProposeState {
  inline ProposeState(std::vector<Eigen::VectorXd> const& biasingPara, std::shared_ptr<ModPiece> const& target, std::shared_ptr<Distribution> const& bias) : biasingPara(biasingPara), target(target), bias(bias) {}

  inline std::shared_ptr<SamplingState> Evaluate(int t) {
    // sample from the biasing distribution
    biasingPara[0] = bias->Sample(std::vector<Eigen::VectorXd>(biasingPara.begin()+1, biasingPara.end()));

    // compute the weight
    const double logweight = target->Evaluate(biasingPara[0])[0](0) - bias->LogDensity(biasingPara);

    // return the sample
    return std::make_shared<SamplingState>(biasingPara[0], std::exp(logweight));
  }

  std::vector<Eigen::VectorXd> biasingPara;

  std::shared_ptr<ModPiece> target;

  std::shared_ptr<Distribution> bias;
};

void ImportanceSampling::RunParallel(std::vector<Eigen::VectorXd>& biasingPara) {
  std::cout << "RUN PARALLEL" << std::endl;
  std::cout << comm->GetSize() << std::endl;

  auto helper = std::make_shared<ProposeState>(biasingPara, target, bias);
  auto proposalQueue = std::make_shared<parcer::Queue<int, std::shared_ptr<SamplingState>, ProposeState> >(helper, comm);

  if( comm->GetRank()==0 ) {
    std::cout << comm->GetRank() << std::endl;
    // submit the work
    std::vector<int> stateIDs(numSamps);

    // submit them to the queue
    for( auto it=stateIDs.begin(); it!=stateIDs.end(); ++it ) {
      int t = 0;
      *it = proposalQueue->SubmitWork(t);
    }

    // retrieve the samples
    /*for( unsigned int i=0; i<numSamps; ++i ) {
      auto state = proposalQueue->GetResult(stateIDs[i]);

      std::cout << state->state[0].transpose() << std::endl;
    }*/
  }
  comm->Barrier();
  //RunSerial(biasingPara);
}
#endif

std::shared_ptr<SampleCollection> ImportanceSampling::RunImpl(std::vector<Eigen::VectorXd> const& x0) {
  // store a copy of the biasing distribution hyper parameters
  std::vector<Eigen::VectorXd> biasingPara = hyperparameters;

  // insert an empty vector for the state
  biasingPara.insert(biasingPara.begin(), Eigen::VectorXd());

#if MUQ_HAS_PARCER
  if( comm->GetSize()==1 ) { RunSerial(biasingPara); }

  RunParallel(biasingPara);
  RunSerial(biasingPara);
#else
  RunSerial(biasingPara);
#endif

  return samples;
}
