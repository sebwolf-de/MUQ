#include "MUQ/SamplingAlgorithms/SLMCMC.h"
#include "MUQ/SamplingAlgorithms/GreedyMLMCMC.h"
#include "MUQ/SamplingAlgorithms/MIMCMC.h"

#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/MHKernel.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/SamplingAlgorithms/CrankNicolsonProposal.h"
#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SubsamplingMIProposal.h"

#include "MUQ/SamplingAlgorithms/MIComponentFactory.h"

#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

#include "MUQ/SamplingAlgorithms/ParallelMIMCMCWorker.h"
#include "MUQ/SamplingAlgorithms/ParallelFixedSamplesMIMCMC.h"

#include "ParallelProblem.h"

#include <ctime>

int main(int argc, char **argv){
  spdlog::set_level(spdlog::level::debug);

  MPI_Init(&argc, &argv);
  auto comm = std::make_shared<parcer::Communicator>();

  // Name trace according to current time stamp
  std::time_t result = std::time(nullptr);
  std::string timestamp = std::asctime(std::localtime(&result));
  auto tracer = std::make_shared<OTF2Tracer>("trace", timestamp);

  pt::ptree pt;
  pt.put("NumSamples_0", 1e3);
  pt.put("NumSamples_1", 5e2);
  pt.put("NumSamples_2", 1e2);
  pt.put("NumSamples_3", 1e2);
  pt.put("MCMC.BurnIn", 100);
  pt.put("MLMCMC.Subsampling", 10);
  pt.put("MLMCMC.Scheduling", true);


  auto componentFactory = std::make_shared<MyMIComponentFactory>(pt);
  StaticLoadBalancingMIMCMC parallelMIMCMC (pt, componentFactory, std::make_shared<RoundRobinStaticLoadBalancer>(), std::make_shared<parcer::Communicator>(), tracer);

  if (comm->GetRank() == 0) {
    parallelMIMCMC.Run();
    Eigen::VectorXd meanQOI = parallelMIMCMC.MeanQOI();
    std::cout << "mean QOI: " << meanQOI.transpose() << std::endl;
  }
  parallelMIMCMC.WriteToFile("parallelMIMCMC.h5");
  parallelMIMCMC.Finalize();
  tracer->write();

  MPI_Finalize();
}
