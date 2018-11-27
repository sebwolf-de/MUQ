#include "MUQ/OptimalExperimentalDesign/Evidence.h"

#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"

#include "MUQ/SamplingAlgorithms/ImportanceSampling.h"
#include "MUQ/SamplingAlgorithms/SampleCollection.h"
#if MUQ_HAS_MPI==1
#include "MUQ/SamplingAlgorithms/DistributedCollection.h"
#endif

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::OptimalExperimentalDesign;

Evidence::Evidence(std::shared_ptr<muq::Modeling::Distribution> const& prior, std::shared_ptr<muq::Modeling::Distribution> const& likelihood, std::shared_ptr<muq::Modeling::Distribution> const& biasing, pt::ptree pt) : Distribution(1, Eigen::VectorXi::Ones(1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), biasing(biasing) {
  CreateGraph(prior, likelihood);
}

Evidence::Evidence(std::shared_ptr<muq::Modeling::Distribution> const& prior, std::shared_ptr<muq::Modeling::Distribution> const& likelihood, std::shared_ptr<muq::Modeling::Distribution> const& biasing, pt::ptree pt, std::shared_ptr<parcer::Communicator> const& comm) : Distribution(1, Eigen::VectorXi::Ones(1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), biasing(biasing), comm(comm) {
  CreateGraph(prior, likelihood);
}

void Evidence::CreateGraph(std::shared_ptr<muq::Modeling::Distribution> const& prior, std::shared_ptr<muq::Modeling::Distribution> const& likelihood) {
  // make a graph to estimate the evidence
  graph = std::make_shared<WorkGraph>();

  // add the input parameters to the graph
  graph->AddNode(std::make_shared<muq::Modeling::IdentityOperator>(1), "design");
  graph->AddNode(std::make_shared<muq::Modeling::IdentityOperator>(1), "data");
  graph->AddNode(std::make_shared<muq::Modeling::IdentityOperator>(1), "parameter");

  // add the densities
  auto priorDens = prior->AsDensity();
  assert(priorDens->numInputs==1);
  graph->AddNode(priorDens, "log prior"); // x
  auto likeDens = likelihood->AsDensity();
  assert(likeDens->numInputs==3);
  graph->AddNode(likeDens, "log likelihood"); // y | x, d
  graph->AddNode(std::make_shared<DensityProduct>(2), "log joint"); // x, y | d

  // connect the prior
  graph->AddEdge("parameter", 0, "log prior", 0);

  // connect the likelihood
  graph->AddEdge("data", 0, "log likelihood", 0);
  graph->AddEdge("parameter", 0, "log likelihood", 1);
  graph->AddEdge("design", 0, "log likelihood", 2);

  // connect the joint
  graph->AddEdge("log prior", 0, "log joint", 0);
  graph->AddEdge("log likelihood", 0, "log joint", 1);
}

double Evidence::LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  // bind the data parameter
  const boost::any& data = inputs[0];
  graph->BindNode("data", std::vector<boost::any>(1, data));

  // bind the design parameter
  const boost::any& design = inputs[1];
  graph->BindNode("design", std::vector<boost::any>(1, design));

  // create the logjoint function
  auto logjoint = graph->CreateModPiece("log joint");

  // create the importance sampler
  pt::ptree pt;
  pt.put("NumSamples", numImportanceSamples);
  auto is = std::make_shared<ImportanceSampling>(logjoint, biasing, pt);

  // compute the importance sampling estimate
#if MUQ_HAS_MPI==1
  std::shared_ptr<SampleCollection> samps;
  if( comm ) {
    auto localSamps = is->Run();
    samps = std::make_shared<DistributedCollection>(localSamps, comm);
  } else {
    samps = is->Run();
  }
#else
  auto samps = is->Run();
#endif
  assert(samps);

  return std::log(samps->Weights().sum()/samps->size());
}
