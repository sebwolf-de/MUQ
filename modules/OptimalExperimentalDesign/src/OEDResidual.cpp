#include "MUQ/OptimalExperimentalDesign/OEDResidual.h"

#include "MUQ/Utilities/RandomGenerator.h"

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
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::Approximation;
using namespace muq::SamplingAlgorithms;
using namespace muq::OptimalExperimentalDesign;

OEDResidual::OEDResidual(std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, pt::ptree pt) : ModPiece(Eigen::Vector2i(likelihood->hyperSizes(0), likelihood->hyperSizes(1)), Eigen::VectorXi::Ones(1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), likelihood(likelihood), evidence(evidence), bruteForce(pt.get<bool>("BruteForce", false)), gamma0(pt.get<double>("InitialThreshold", 1.0)), radius0(pt.get<double>("InitialRadius", 1.0)) {
  CreateGraph();

  if( !bruteForce ) { pt_regression = pt.get_child(pt.get<std::string>("LocalRegression")); }
}

#if MUQ_HAS_PARCER==1
OEDResidual::OEDResidual(std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, pt::ptree pt, std::shared_ptr<parcer::Communicator> const& comm) : ModPiece(Eigen::Vector2i(likelihood->varSize, likelihood->hyperSizes(1)), Eigen::VectorXi::Ones(1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), likelihood(likelihood), evidence(evidence), bruteForce(pt.get<bool>("BruteForce", false)), gamma0(pt.get<double>("InitialThreshold", 1.0)), radius0(pt.get<double>("InitialRadius", 1.0)), comm(comm) {
  CreateGraph();

  if( !bruteForce ) { pt_regression = pt.get_child(pt.get<std::string>("LocalRegression")); }
}
#endif

void OEDResidual::CreateGraph() {
  // make a graph to estimate the evidence
  graph = std::make_shared<WorkGraph>();

  // add the input parameters to the graph
  graph->AddNode(std::make_shared<IdentityOperator>(inputSizes(1)), "design");
  graph->AddNode(std::make_shared<IdentityOperator>(likelihood->varSize), "data");
  graph->AddNode(std::make_shared<IdentityOperator>(inputSizes(0)), "parameter");

  // add the densities
  graph->AddNode(likelihood->AsDensity(), "log likelihood"); // y | x, d

  // the log difference between the evidence and the likelihood
  graph->AddNode(std::make_shared<LogDifference>(likelihood, evidence), "log difference");

  // connect the likelihood
  graph->AddEdge("data", 0, "log likelihood", 0);
  graph->AddEdge("parameter", 0, "log likelihood", 1);
  graph->AddEdge("design", 0, "log likelihood", 2);

  // connect the difference
  graph->AddEdge("data", 0, "log difference", 0);
  graph->AddEdge("parameter", 0, "log difference", 1);
  graph->AddEdge("design", 0, "log difference", 2);
}

void OEDResidual::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) {
  const boost::any& parameter = inputs[0];
  graph->BindNode("parameter", std::vector<boost::any>(1, parameter));

  const boost::any& design = inputs[1];
  graph->BindNode("design", std::vector<boost::any>(1, design));

  auto loglikelihood = graph->CreateModPiece("log likelihood");

  pt::ptree pt;
  pt.put("NumSamples", numImportanceSamples);
  auto is = std::make_shared<ImportanceSampling>(likelihood, std::vector<Eigen::VectorXd>({inputs[0], inputs[1]}), pt);

std::shared_ptr<SampleCollection> samps, localSamps;
#if MUQ_HAS_MPI==1
  if( comm ) {
    localSamps = is->Run();
    samps = std::make_shared<DistributedCollection>(localSamps, comm);
  } else {
    samps = is->Run();
    localSamps = samps;
  }
#else
  samps = is->Run();
  localSamps = samps;
#endif
  assert(samps);

  auto diff = graph->CreateModPiece("log difference");
  assert(diff);

  outputs.resize(1);
  if( bruteForce ) {
    outputs[0] = samps->ExpectedValue(diff);
    return;
  }
  
#if MUQ_HAS_MPI==1
  auto reg = comm ? std::make_shared<LocalRegression>(diff, pt_regression, comm) : std::make_shared<LocalRegression>(diff, pt_regression);
#else
  auto reg = std::make_shared<LocalRegression>(diff, pt_regression);
#endif
  const double threshold = std::log(gamma0)-0.5*std::log((double)samps->size());
  totalRefinements = 0;

  for( unsigned int i=0; i<localSamps->size(); ++i ) {
    const Eigen::VectorXd& xd = localSamps->at(i)->state[0];

    // make sure there are enough points in the cache
    while( reg->CacheSize()<reg->kn ) {
      RandomlyRefineNear(reg, xd, radius0);
    }

    const std::pair<double, double> error = reg->ErrorIndicator(xd);

    //std::cout << "OED RESIDUAL, step: " << i << " error: " << error.first << " threshold: " << threshold << std::endl;
    //if( error.first>threshold ) {
      //std::cout << "TRIGGERED!" << std::endl;
      const std::tuple<Eigen::VectorXd, double, unsigned int> lambda = reg->PoisednessConstant(xd);
        RefineAt(reg, std::get<0>(lambda), error.second);
      //}
  }

  std::cout << "OED RESIDUAL, total refinements: " << totalRefinements << std::endl;
  outputs[0] = samps->ExpectedValue(reg);
}

void OEDResidual::RefineAt(std::shared_ptr<LocalRegression> const& reg, Eigen::VectorXd const& pnt, double const radius) {
  assert(reg);
  if( reg->InCache(pnt) ) { return RandomlyRefineNear(reg, pnt, radius); }

  reg->Add(pnt);
  ++totalRefinements;
}


void OEDResidual::RandomlyRefineNear(std::shared_ptr<LocalRegression> const& reg, Eigen::VectorXd const& xd, double const radius) {
  assert(reg);

  while( true ) {
    // choose a random point
    Eigen::VectorXd point = RandomGenerator::GetNormal(xd.size());
    point *= RandomGenerator::GetUniform()*radius/point.norm();
    point += xd;
    if( !reg->InCache(point) ) {
      reg->Add(point);
      ++totalRefinements;
      break;
    }
  }
}
