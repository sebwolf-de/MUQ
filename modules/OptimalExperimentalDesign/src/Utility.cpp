#include "MUQ/OptimalExperimentalDesign/Utility.h"

#include "MUQ/Modeling/ScaleVector.h"
#include "MUQ/Modeling/SplitVector.h"
#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"

#include "MUQ/Optimization/ModPieceCostFunction.h"

#include "MUQ/SamplingAlgorithms/ImportanceSampling.h"
#include "MUQ/SamplingAlgorithms/SampleCollection.h"
#if MUQ_HAS_MPI==1
#include "MUQ/SamplingAlgorithms/DistributedCollection.h"
#endif
#include "MUQ/OptimalExperimentalDesign/LogDifference.h"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Optimization;
using namespace muq::SamplingAlgorithms;
using namespace muq::OptimalExperimentalDesign;

Utility::Utility(std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, std::shared_ptr<Distribution> const& biasing, boost::property_tree::ptree pt) : ModPiece(Eigen::VectorXi::Ones(1), Eigen::VectorXi::Ones(1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), evidence(evidence), biasing(biasing) {
  CreateGraph(prior, likelihood, evidence);
}

#if MUQ_HAS_PARCER==1
Utility::Utility(std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, std::shared_ptr<Distribution> const& biasing, boost::property_tree::ptree pt, std::shared_ptr<parcer::Communicator> const& comm) : ModPiece(Eigen::VectorXi::Ones(1), Eigen::VectorXi::Ones(1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), evidence(evidence), biasing(biasing), comm(comm) {
  CreateGraph(prior, likelihood, evidence);
}
#endif

void Utility::CreateGraph(std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence) {
  // make a graph to estimate the evidence
  graph = std::make_shared<WorkGraph>();

  // add the input parameters to the graph
  graph->AddNode(std::make_shared<IdentityOperator>(1), "design");
  graph->AddNode(std::make_shared<IdentityOperator>(1), "data");
  graph->AddNode(std::make_shared<IdentityOperator>(1), "parameter");

  // add a node to split the parameters
  graph->AddNode(std::make_shared<SplitVector>(Eigen::Vector2i(0, 1), Eigen::VectorXi::Ones(2), 2), "parameter-data");

  // add the densities
  graph->AddNode(prior->AsDensity(), "log prior"); // x
  graph->AddNode(likelihood->AsDensity(), "log likelihood"); // y | x, d
  graph->AddNode(std::make_shared<DensityProduct>(2), "log joint"); // x, y | d
  graph->AddNode(std::make_shared<ScaleVector>(-1.0, 1), "negative log joint");

  // the log difference between the evidence and the likelihood
  graph->AddNode(std::make_shared<LogDifference>(likelihood, evidence), "log difference");

  // connect the parameters
  graph->AddEdge("parameter-data", 0, "parameter", 0);
  graph->AddEdge("parameter-data", 1, "data", 0);

  // connect the prior
  graph->AddEdge("parameter", 0, "log prior", 0);

  // connect the likelihood
  graph->AddEdge("data", 0, "log likelihood", 0);
  graph->AddEdge("parameter", 0, "log likelihood", 1);
  graph->AddEdge("design", 0, "log likelihood", 2);

  // connect the joint
  graph->AddEdge("log prior", 0, "log joint", 0);
  graph->AddEdge("log likelihood", 0, "log joint", 1);
  graph->AddEdge("log joint", 0, "negative log joint", 0);

  // connect the difference
  graph->AddEdge("parameter", 0, "log difference", 0);
  graph->AddEdge("data", 0, "log difference", 1);
  graph->AddEdge("design", 0, "log difference", 2);
}

void Utility::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  const boost::any& design = inputs[0];
  graph->BindNode("design", std::vector<boost::any>(1, design));

  auto logjoint = std::make_shared<ModPieceCostFunction>(graph->CreateModPiece("log joint"));
  assert(logjoint);
  assert(biasing);

  pt::ptree pt;
  pt.put<unsigned int>("ImportanceSampling.NumSamples", numImportanceSamples);
  auto is = std::make_shared<ImportanceSampling>(logjoint, biasing, std::vector<Eigen::VectorXd>(1, inputs[0]), pt.get_child("ImportanceSampling"));

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

  auto diff = graph->CreateModPiece("log difference");
  assert(diff);

  outputs.resize(1);
  outputs[0] = samps->ExpectedValue(diff);
}
