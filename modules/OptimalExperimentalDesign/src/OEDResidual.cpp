#include "MUQ/OptimalExperimentalDesign/OEDResidual.h"

#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"
#include "MUQ/Modeling/SplitVector.h"
#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/Distributions/Density.h"

#include "MUQ/SamplingAlgorithms/ImportanceSampling.h"
#include "MUQ/SamplingAlgorithms/SampleCollection.h"
#if MUQ_HAS_MPI==1
#include "MUQ/SamplingAlgorithms/DistributedCollection.h"
#endif

#include "MUQ/OptimalExperimentalDesign/LogDifference.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::OptimalExperimentalDesign;

OEDResidual::OEDResidual(std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, pt::ptree pt) : ModPiece(Eigen::Vector2i(likelihood->hyperSizes(0), likelihood->hyperSizes(1)), Eigen::VectorXi::Ones(1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), evidence(evidence), biasing(likelihood) {
  CreateGraph(likelihood, evidence);
}

OEDResidual::OEDResidual(std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, std::shared_ptr<Distribution> const& biasing, pt::ptree pt) : ModPiece(Eigen::Vector2i(likelihood->hyperSizes(0), likelihood->hyperSizes(1)), Eigen::VectorXi::Ones(1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), evidence(evidence), biasing(biasing) {
  CreateGraph(likelihood, evidence);
}

#if MUQ_HAS_PARCER==1
OEDResidual::OEDResidual(std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, pt::ptree pt, std::shared_ptr<parcer::Communicator> const& comm) : ModPiece(Eigen::Vector2i(likelihood->varSize, likelihood->hyperSizes(1)), Eigen::VectorXi::Ones(1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), evidence(likelihood), biasing(biasing), comm(comm) {
  CreateGraph(likelihood, evidence);
}

OEDResidual::OEDResidual(std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, std::shared_ptr<Distribution> const& biasing, pt::ptree pt, std::shared_ptr<parcer::Communicator> const& comm) : ModPiece(Eigen::Vector2i(likelihood->varSize, likelihood->hyperSizes(1)), Eigen::VectorXi::Ones(1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), evidence(evidence), biasing(biasing), comm(comm) {
  CreateGraph(likelihood, evidence);
}
#endif

void OEDResidual::CreateGraph(std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence) {
  // make a graph to estimate the evidence
  graph = std::make_shared<WorkGraph>();

  // add the input parameters to the graph
  graph->AddNode(std::make_shared<IdentityOperator>(inputSizes(1)), "design");
  graph->AddNode(std::make_shared<IdentityOperator>(likelihood->varSize), "data");
  graph->AddNode(std::make_shared<IdentityOperator>(inputSizes(0)), "parameter");

  // add the densities
  graph->AddNode(likelihood->AsDensity(), "log likelihood"); // y | x, d

  // the log difference between the evidence and the likelihood
  graph->AddNode(std::make_shared<LogDifference>(evidence), "log difference");

  // connect the likelihood
  graph->AddEdge("data", 0, "log likelihood", 0);
  graph->AddEdge("parameter", 0, "log likelihood", 1);
  graph->AddEdge("design", 0, "log likelihood", 2);

  // connect the difference
  graph->AddEdge("data", 0, "log difference", 0);
  graph->AddEdge("design", 0, "log difference", 1);
}

void OEDResidual::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) {
  const boost::any& parameter = inputs[0];
  graph->BindNode("parameter", std::vector<boost::any>(1, parameter));

  const boost::any& design = inputs[1];
  graph->BindNode("design", std::vector<boost::any>(1, design));

  auto loglikelihood = graph->CreateModPiece("log likelihood");

  pt::ptree pt;
  pt.put("NumSamples", numImportanceSamples);
  auto is = std::make_shared<ImportanceSampling>(loglikelihood, biasing, std::vector<Eigen::VectorXd>({inputs[0], inputs[1]}), pt);

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
  
  auto diff = graph->CreateModPiece("log difference", std::vector<std::string>({"data", "log difference"}));
  assert(diff);

  outputs.resize(1);
  outputs[0] = samps->ExpectedValue(diff, std::vector<std::string>({"log target"}));
}
