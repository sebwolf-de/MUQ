#include "MUQ/OptimalExperimentalDesign/Utility.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/Modeling/SplitVector.h"
#include "MUQ/Modeling/CombineVectors.h"
#include "MUQ/Modeling/ConstantVector.h"
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
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::Optimization;
using namespace muq::Approximation;
using namespace muq::SamplingAlgorithms;
using namespace muq::OptimalExperimentalDesign;

Utility::Utility(std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, pt::ptree pt) : ModPiece(Eigen::VectorXi::Constant(1, likelihood->hyperSizes(1)), Eigen::VectorXi::Constant(pt.get<bool>("RunningEstimate", false) ? pt.get<unsigned int>("NumImportanceSamples") : 1, 1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), biasing(prior), runningEstimate(pt.get<bool>("RunningEstimate", false)) {
  CreateGraph(prior, likelihood, evidence);
}

Utility::Utility(std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, std::shared_ptr<Distribution> const& biasing, pt::ptree pt) : ModPiece(Eigen::VectorXi::Constant(1, likelihood->hyperSizes(1)), Eigen::VectorXi::Constant(pt.get<bool>("RunningEstimate", false) ? pt.get<unsigned int>("NumImportanceSamples") : 1, 1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), biasing(biasing), runningEstimate(pt.get<bool>("RunningEstimate", false)) {
  CreateGraph(prior, likelihood, evidence);
}

Utility::Utility(std::shared_ptr<Distribution> const& prior, std::shared_ptr<OEDResidual> const& resid, pt::ptree pt) : ModPiece(Eigen::VectorXi::Constant(1, resid->inputSizes(1)), Eigen::VectorXi::Constant(pt.get<bool>("RunningEstimate", false) ? pt.get<unsigned int>("NumImportanceSamples") : 1, 1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), gamma0(pt.get<double>("InitialThreshold")), radius0(pt.get<double>("InitialRadius")), biasing(prior), runningEstimate(pt.get<bool>("RunningEstimate", false)) {
  CreateGraph(prior, resid, pt);
}

Utility::Utility(std::shared_ptr<muq::Modeling::Distribution> const& prior, std::shared_ptr<OEDResidual> const& resid, std::shared_ptr<muq::Modeling::Distribution> const& biasing, pt::ptree pt) : ModPiece(Eigen::VectorXi::Constant(1, resid->inputSizes(1)), Eigen::VectorXi::Constant(pt.get<bool>("RunningEstimate", false) ? pt.get<unsigned int>("NumImportanceSamples") : 1, 1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), gamma0(pt.get<double>("InitialThreshold")), radius0(pt.get<double>("InitialRadius")), biasing(biasing), runningEstimate(pt.get<bool>("RunningEstimate", false)) {
  CreateGraph(prior, resid, pt);
}

#if MUQ_HAS_PARCER==1
Utility::Utility(std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, pt::ptree pt, std::shared_ptr<parcer::Communicator> const& comm) : ModPiece(Eigen::VectorXi::Constant(1, likelihood->hyperSizes(1)), Eigen::VectorXi::Constant(pt.get<bool>("RunningEstimate", false) ? pt.get<unsigned int>("NumImportanceSamples") : 1, 1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), biasing(prior), runningEstimate(pt.get<bool>("RunningEstimate", false)), comm(comm) {
  if( comm->GetSize()>1 ) { assert(!runningEstimate); }
  CreateGraph(prior, likelihood, evidence);
}

Utility::Utility(std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence, std::shared_ptr<Distribution> const& biasing, pt::ptree pt, std::shared_ptr<parcer::Communicator> const& comm) : ModPiece(Eigen::VectorXi::Constant(1, likelihood->hyperSizes(1)), Eigen::VectorXi::Constant(pt.get<bool>("RunningEstimate", false) ? pt.get<unsigned int>("NumImportanceSamples") : 1, 1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), biasing(biasing), runningEstimate(pt.get<bool>("RunningEstimate", false)), comm(comm) {
  if( comm->GetSize()>1 ) { assert(!runningEstimate); }
  CreateGraph(prior, likelihood, evidence);
}

Utility::Utility(std::shared_ptr<Distribution> const& prior, std::shared_ptr<OEDResidual> const& resid, pt::ptree pt, std::shared_ptr<parcer::Communicator> const& comm) : ModPiece(Eigen::VectorXi::Constant(1, resid->inputSizes(1)), Eigen::VectorXi::Constant(pt.get<bool>("RunningEstimate", false) ? pt.get<unsigned int>("NumImportanceSamples") : 1, 1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), gamma0(pt.get<double>("InitialThreshold")), radius0(pt.get<double>("InitialRadius")), biasing(prior), runningEstimate(pt.get<bool>("RunningEstimate", false)), comm(comm) {
  if( comm->GetSize()>1 ) { assert(!runningEstimate); }
  CreateGraph(prior, resid, pt);
}

Utility::Utility(std::shared_ptr<Distribution> const& prior, std::shared_ptr<OEDResidual> const& resid, std::shared_ptr<Distribution> const& biasing, pt::ptree pt, std::shared_ptr<parcer::Communicator> const& comm) : ModPiece(Eigen::VectorXi::Constant(1, resid->inputSizes(1)), Eigen::VectorXi::Constant(pt.get<bool>("RunningEstimate", false) ? pt.get<unsigned int>("NumImportanceSamples") : 1, 1)), numImportanceSamples(pt.get<unsigned int>("NumImportanceSamples")), gamma0(pt.get<double>("InitialThreshold")), radius0(pt.get<double>("InitialRadius")), biasing(biasing), runningEstimate(pt.get<bool>("RunningEstimate", false)), comm(comm) {
  if( comm->GetSize()>1 ) { assert(!runningEstimate); }
  CreateGraph(prior, resid, pt);
}
#endif

void Utility::CreateGraph(std::shared_ptr<muq::Modeling::Distribution> const& prior, std::shared_ptr<OEDResidual> const& resid, pt::ptree pt) {
  // make a graph to estimate the evidence
  graph = std::make_shared<WorkGraph>();

  // add the input parameters to the graph
  graph->AddNode(std::make_shared<IdentityOperator>(resid->inputSizes(1)), "design");
  graph->AddNode(std::make_shared<IdentityOperator>(resid->inputSizes(0)), "parameter");

  // add the densities
  graph->AddNode(prior->AsDensity(), "log prior"); // x

  // add a node to split the parameters
  graph->AddNode(std::make_shared<SplitVector>(Eigen::Vector2i(0, resid->inputSizes(0)), Eigen::Vector2i(resid->inputSizes(0), resid->inputSizes(1)), resid->inputSizes(0)+resid->inputSizes(1)), "parameter-design");

  // add the residual
  graph->AddNode(resid, "residual");

  // link the parameters and design variables
  graph->AddEdge("parameter-design", 0, "parameter", 0);
  graph->AddEdge("parameter-design", 1, "design", 0);

  // link the residual
  graph->AddEdge("parameter", 0, "residual", 0);
  graph->AddEdge("design", 0, "residual", 1);

  // link the prior
  graph->AddEdge("parameter", 0, "log prior", 0);

  // create the local regression
#if MUQ_HAS_MPI==1
  if( comm ) {
    reg = std::make_shared<LocalRegression>(graph->CreateModPiece("residual"), pt.get_child(pt.get<std::string>("LocalRegression")), comm);
  } else {
    reg = std::make_shared<LocalRegression>(graph->CreateModPiece("residual"), pt.get_child(pt.get<std::string>("LocalRegression")));
  }
#else
  reg = std::make_shared<LocalRegression>(graph->CreateModPiece("residual"), pt.get_child(pt.get<std::string>("LocalRegression")));
#endif

  // remove the combined parameter/design node
  graph->RemoveNode("parameter-design");
}

void Utility::CreateGraph(std::shared_ptr<Distribution> const& prior, std::shared_ptr<Distribution> const& likelihood, std::shared_ptr<Distribution> const& evidence) {
  // make a graph to estimate the evidence
  graph = std::make_shared<WorkGraph>();

  // add the input parameters to the graph
  graph->AddNode(std::make_shared<IdentityOperator>(likelihood->hyperSizes(1)), "design");
  graph->AddNode(std::make_shared<IdentityOperator>(likelihood->varSize), "data");
  graph->AddNode(std::make_shared<IdentityOperator>(likelihood->hyperSizes(0)), "parameter");

  // add a node to split the parameters
  graph->AddNode(std::make_shared<SplitVector>(Eigen::Vector2i(0, likelihood->hyperSizes(0)), Eigen::Vector2i(likelihood->hyperSizes(0), likelihood->varSize), likelihood->hyperSizes(0)+likelihood->varSize), "parameter-data");

  // add the densities
  graph->AddNode(prior->AsDensity(), "log prior"); // x
  graph->AddNode(likelihood->AsDensity(), "log likelihood"); // y | x, d
  graph->AddNode(std::make_shared<DensityProduct>(2), "log joint"); // x, y | d

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

  // connect the difference
  graph->AddEdge("data", 0, "log difference", 0);
  graph->AddEdge("parameter", 0, "log difference", 1);
  graph->AddEdge("design", 0, "log difference", 2);
}

void Utility::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  if( reg ) {
    EvaluateSurrogate(inputs);
    return;
  }

  EvaluateBruteForce(inputs);
}

void Utility::RandomlyRefineNear(Eigen::VectorXd const& xd, double const radius) {
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

void Utility::RefineAt(Eigen::VectorXd const& pnt, double const radius) {
  assert(reg);
  if( reg->InCache(pnt) ) { return RandomlyRefineNear(pnt, radius); }

  reg->Add(pnt);
  ++totalRefinements;
}

void Utility::EvaluateSurrogate(ref_vector<Eigen::VectorXd> const& inputs) {
  //std::cout << "EVALUATE SURROGATE" << std::endl;
  auto logprior = graph->CreateModPiece("log prior");
  assert(logprior);
  assert(biasing);

  // run the importance sampler
  pt::ptree pt;
  pt.put("NumSamples", numImportanceSamples);
  auto is = std::make_shared<ImportanceSampling>(logprior, biasing, std::vector<Eigen::VectorXd>(1, inputs[0]), pt);

  #if MUQ_HAS_MPI==1
    std::shared_ptr<SampleCollection> samps;
    std::shared_ptr<SampleCollection> localSamps;
    if( comm ) {
      localSamps = is->Run();
      samps = std::make_shared<DistributedCollection>(localSamps, comm);
    } else {
      samps = is->Run();
      localSamps = samps;
    }
  #else
    auto samps = is->Run();
    auto localSamps = samps;
  #endif
  assert(samps);
  assert(localSamps);

  const unsigned int n2 = samps->size();
  //const double threshold = gamma0*std::pow((double)n2, -0.5);
  const double threshold = std::log(gamma0)-0.5*std::log((double)n2);

  // loop through the samples
  Eigen::VectorXd xd(inputSizes(0)+logprior->inputSizes(0));
  xd.tail(inputSizes(0)) = inputs[0].get();

  /*auto g = std::make_shared<WorkGraph>();

  g->AddNode(reg, "residual approx");
  g->AddNode(std::make_shared<CombineVectors>(Eigen::Vector2i(logprior->inputSizes(0), inputSizes(0))), "combine");
  g->AddNode(std::make_shared<IdentityOperator>(logprior->inputSizes(0)), "parameter");
  g->AddNode(std::make_shared<ConstantVector>(std::vector<Eigen::VectorXd>(1, inputs[0])), "design");

  g->AddEdge("parameter", 0, "combine", 0);
  g->AddEdge("design", 0, "combine", 1);
  g->AddEdge("combine", 0, "residual approx", 0);

  auto mod = g->CreateModPiece("residual approx");*/

  if( runningEstimate ) { outputs.resize(localSamps->size()); }

  totalRefinements = 0;
  runningRefinements = Eigen::VectorXi::Zero(localSamps->size());
  for( unsigned int i=0; i<localSamps->size(); ++i ) {
    // the current point for the surrogate model
    xd.head(logprior->inputSizes(0)) = localSamps->at(i)->state[0];

    // make sure there are enough points in the cache
    while( reg->CacheSize()<reg->kn ) {
      //if( comm )
        //std::cout << "process " << comm->GetRank() << " is refining, needs more neighbors: step " << i << std::endl;
      //else
      //std::cout << "refining, needs more neighbors; step " << i << std::endl;
      RandomlyRefineNear(xd, radius0);
    }

    const std::pair<double, double> error = reg->ErrorIndicator(xd);
    //std::cout << "thershold: " << std::exp(threshold) << " error: " << std::exp(error.first) << "  radius: " << error.second << std::endl;
    //std::cout << "total refinements: " << totalRefinements << std::endl;
    //const std::tuple<Eigen::VectorXd, double, unsigned int> lambda = reg->PoisednessConstant(xd);
    //std::cout << "\tpoisedness constant: " << std::get<1>(lambda) << " step: " << i << std::endl;

    /*if( std::get<1>(lambda)>100.0 ) {
      std::cout << "\tprocess " << comm->GetRank() << " poisedness constant: " << std::get<1>(lambda) << " step: " << i << std::endl;
    }*/

    if( error.first>threshold /*|| std::get<1>(lambda)>lammax*/ ) {
      const std::tuple<Eigen::VectorXd, double, unsigned int> lambda = reg->PoisednessConstant(xd);
      //std::cout << "\tpoisedness constant: " << std::get<1>(lambda) << " step: " << i << std::endl;

      /*if( comm ) {
        std::cout << "process " << comm->GetRank() << " is refining, indicator: " << error.first << " threshold: " << threshold << " step: " << i << std::endl;
        std::cout << "\tpoisedness: " << std::get<1>(lambda) << std::endl;
      } else {
        std::cout << "refining, indicator: " << error.first << " threshold: " << threshold << " step: " << i << std::endl;
      }*/
      RefineAt(std::get<0>(lambda), error.second);

      //error = reg->ErrorIndicator(xd);
      //std::cout << "\tthershold: " << threshold << " error: " << error.first << "  radius: " << error.second << std::endl;
    }

    /*if( runningEstimate ) {
      auto g = std::make_shared<WorkGraph>();

      g->AddNode(reg, "residual approx");
      g->AddNode(std::make_shared<CombineVectors>(Eigen::Vector2i(logprior->inputSizes(0), inputSizes(0))), "combine");
      g->AddNode(std::make_shared<IdentityOperator>(logprior->inputSizes(0)), "parameter");
      g->AddNode(std::make_shared<ConstantVector>(std::vector<Eigen::VectorXd>(1, inputs[0])), "design");

      g->AddEdge("parameter", 0, "combine", 0);
      g->AddEdge("design", 0, "combine", 1);
      g->AddEdge("combine", 0, "residual approx", 0);

      auto mod = g->CreateModPiece("residual approx");

      double totWeight = 0.0;
      Eigen::VectorXd ubreve = Eigen::VectorXd::Zero(1);
      for( unsigned int j=0; j<i; ++j ) {
        auto s = localSamps->at(j);
        totWeight += s->weight;
        if( i==0 ) {
          ubreve = mod->Evaluate(s->state) [0];
        } else {
          ubreve = ((totWeight-s->weight)*ubreve + s->weight*mod->Evaluate(s->state) [0])/totWeight;
        }
      }
      std::cout << "ubreve: " << ubreve << std::endl;
      outputs[i] = ubreve;
    }*/

    runningRefinements(i) = totalRefinements;
  }

//  if( !runningEstimate ) {
    auto g = std::make_shared<WorkGraph>();

    g->AddNode(reg, "residual approx");
    g->AddNode(std::make_shared<CombineVectors>(Eigen::Vector2i(logprior->inputSizes(0), inputSizes(0))), "combine");
    g->AddNode(std::make_shared<IdentityOperator>(logprior->inputSizes(0)), "parameter");
    g->AddNode(std::make_shared<ConstantVector>(std::vector<Eigen::VectorXd>(1, inputs[0])), "design");

    g->AddEdge("parameter", 0, "combine", 0);
    g->AddEdge("design", 0, "combine", 1);
    g->AddEdge("combine", 0, "residual approx", 0);

    auto mod = g->CreateModPiece("residual approx");

    outputs.resize(1);
    outputs[0] = samps->ExpectedValue(mod);
  //}

  /*Eigen::VectorXd d = inputs[0].get();
  Eigen::VectorXd expected = Eigen::VectorXd::Zero(1);
  for( unsigned int i=0; i<samps->size(); ++i ) {
    std::cout << i+1 << " of " << samps->size() << std::endl;
    // the current point for the surrogate model
    xd.head(logprior->inputSizes(0)) = samps->at(i)->state[0];
    Eigen::VectorXd x = samps->at(i)->state[0];

    //const Eigen::VectorXd truth = TEST->Evaluate(xd) [0];
    const Eigen::VectorXd truth = TEST->Evaluate(x, d) [0];
    expected += truth;

    const Eigen::VectorXd approx = reg->Evaluate(xd) [0];
    std::cout << "truth: " << truth.transpose() << std::endl;
    std::cout << "approx: " << approx.transpose() << std::endl;
    std::cout << std::endl;
  }*/
  //outputs.resize(1);
  //outputs[0] = (Eigen::VectorXd)(expected/samps->size());
}

void Utility::EvaluateBruteForce(ref_vector<Eigen::VectorXd> const& inputs) {
  // bind the design node
  const boost::any& design = inputs[0];
  graph->BindNode("design", std::vector<boost::any>(1, design));

  auto logjoint = graph->CreateModPiece("log joint");
  assert(logjoint);
  assert(biasing);

  pt::ptree pt;
  pt.put("NumSamples", numImportanceSamples);
  auto is = std::make_shared<ImportanceSampling>(logjoint, biasing, std::vector<Eigen::VectorXd>(1, inputs[0]), pt);

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

  auto diff = graph->CreateModPiece("log difference", std::vector<std::string>({"parameter-data"}));
  assert(diff);

  if( runningEstimate ) {
    outputs = samps->RunningExpectedValue(diff);
  } else {
    outputs.resize(1);
    outputs[0] = samps->ExpectedValue(diff);
  }
}

unsigned int Utility::TotalRefinements() const { return totalRefinements; }

Eigen::VectorXi Utility::RunningRefinements() const { return runningRefinements; }
