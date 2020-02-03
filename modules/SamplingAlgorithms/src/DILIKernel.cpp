#include "MUQ/SamplingAlgorithms/DILIKernel.h"

#include "MUQ/Modeling/LinearAlgebra/HessianOperator.h"
#include "MUQ/Modeling/LinearAlgebra/GaussNewtonOperator.h"
#include "MUQ/Modeling/LinearAlgebra/GaussianOperator.h"
#include "MUQ/Modeling/LinearAlgebra/LOBPCG.h"
#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"
#include "MUQ/Modeling/SumPiece.h"
#include "MUQ/Modeling/ModGraphPiece.h"

using namespace muq::SamplingAlgorithms;
using namespace muq::Modeling;


Eigen::MatrixXd CSProjector::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
  return x - (*U)*W->transpose() * x;
}

Eigen::MatrixXd CSProjector::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
  return x - (*W)*U->transpose() * x;
}

Eigen::MatrixXd LIS2Full::Apply(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
  return (*U) * L->triangularView<Eigen::Lower>() * x;
}

Eigen::MatrixXd LIS2Full::ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x)
{
  return L->triangularView<Eigen::Lower>().transpose() * U->transpose() * x;
}



DILIKernel::DILIKernel(boost::property_tree::ptree       const& pt,
                       std::shared_ptr<AbstractSamplingProblem> problem) : TransitionKernel(pt,problem)
{
 assert(false);
}


DILIKernel::DILIKernel(boost::property_tree::ptree                  const& pt,
                       std::shared_ptr<AbstractSamplingProblem>            problem,
                       std::shared_ptr<muq::Modeling::ModPiece>     const& likelihoodIn,
                       std::shared_ptr<muq::Modeling::GaussianBase> const& priorIn) : TransitionKernel(pt,problem),
                                                                                      lisKernelOpts(pt.get_child(pt.get<std::string>("LIS Block"))),
                                                                                      csKernelOpts(pt.get_child(pt.get<std::string>("CS Block"))),
                                                                                      logLikelihood(likelihoodIn),
                                                                                      prior(priorIn),
                                                                                      hessType(pt.get("HessianType","Exact"))
{
  try{
    std::string blockName = pt.get<std::string>("Eigensolver Block");
    lobpcgOpts = pt.get_child(blockName);
  }catch(boost::property_tree::ptree_bad_path){
    // Do nothing, just leave the solver options ptree empty
  }
}

void DILIKernel::PostStep(unsigned int const t,
                          std::vector<std::shared_ptr<SamplingState>> const& state)
{
}

std::vector<std::shared_ptr<SamplingState>> DILIKernel::Step(unsigned int const t,
                                                             std::shared_ptr<SamplingState> prevState)
{
  if(lisU==nullptr){
    CreateLIS(prevState->state, blockInd, lobpcgOpts);
  }

  std::vector<Eigen::VectorXd> splitVec(2);
  splitVec.at(0) = lisL->triangularView<Eigen::Lower>().solve(lisW->transpose()*prevState->state.at(blockInd));
  splitVec.at(1) = prevState->state.at(blockInd);

  std::shared_ptr<SamplingState> splitState = std::make_shared<SamplingState>(splitVec);
  splitState->meta = prevState->meta;

  // Take the alternating Metropolis-in-Gibbs steps for the LIS and CS
  std::vector<std::shared_ptr<SamplingState>> nextSteps = lisKernel->Step(t, splitState);
  std::vector<std::shared_ptr<SamplingState>> nextSteps2 = csKernel->Step(t+nextSteps.size(), nextSteps.at(nextSteps.size()-1));

  // Combine the states and return the results
  std::vector<std::shared_ptr<SamplingState>> output(nextSteps.size() + nextSteps2.size());
  Eigen::VectorXd x;
  for(unsigned int i=0; i<nextSteps.size(); ++i){
    x = (*lisU)*lisL->triangularView<Eigen::Lower>() * nextSteps.at(i)->state.at(0);
    x += nextSteps.at(i)->state.at(1) - (*lisU) * lisW->transpose() * nextSteps.at(i)->state.at(1);
    output.at(i) = std::make_shared<SamplingState>(prevState->state);
    output.at(i)->state.at(blockInd) = x;
  }
  for(unsigned int i=0; i<nextSteps2.size(); ++i){
    x = (*lisU)*lisL->triangularView<Eigen::Lower>() * nextSteps2.at(i)->state.at(0);
    x += nextSteps2.at(i)->state.at(1) - (*lisU) * lisW->transpose() * nextSteps2.at(i)->state.at(1);
    output.at(i+nextSteps.size()) = std::make_shared<SamplingState>(prevState->state);
    output.at(i+nextSteps.size())->state.at(blockInd) = x;
  }

  return output;
}

void DILIKernel::PrintStatus(std::string prefix) const
{
  std::string newPrefix = prefix + " LIS: ";
  lisKernel->PrintStatus(newPrefix);
  newPrefix = prefix + " CS: ";
  csKernel->PrintStatus(newPrefix);
}

std::shared_ptr<muq::Modeling::ModPiece> DILIKernel::ExtractLikelihood(std::shared_ptr<muq::Modeling::ModPiece> const& tgtGraphPiece,
                                                                       std::string                              const& nodeName)
{
assert(false);
}


std::shared_ptr<muq::Modeling::GaussianBase> DILIKernel::ExtractPrior(std::shared_ptr<muq::Modeling::ModPiece> const& tgtGraphPiece,
                                                                      std::string                              const& nodeName)
{
assert(false);
}

void DILIKernel::CreateLIS(std::vector<Eigen::VectorXd> const& currState,
                           unsigned int                        blockId,
                           boost::property_tree::ptree  const& solverOpts)
{
  // First, set up the hessian operator for the log-likelihood
  std::shared_ptr<LinearOperator> hessOp;
  if(hessType=="Exact"){
    hessOp = std::make_shared<HessianOperator>(logLikelihood, currState, 0, blockId, blockId, Eigen::VectorXd::Ones(1), -1.0, 0.0);
  }else if(hessType=="GaussNewton"){
    hessOp = std::make_shared<GaussNewtonOperator>(forwardModel, noiseModel, currState, blockId);
  }else{
    std::cerr << "\nERROR: Unrecognized Hessian type.  Options are \"Exact\" or \"GaussNewton\".\n\n";
  }

  // Set up the prior precision operator
  std::shared_ptr<LinearOperator> precOp = std::make_shared<GaussianOperator>(prior, Gaussian::Precision);
  std::shared_ptr<LinearOperator> covOp = std::make_shared<GaussianOperator>(prior, Gaussian::Covariance);

  // Solve the generalized Eigenvalue problem using LOBPCG
  LOBPCG solver(solverOpts);
  solver.compute(hessOp, precOp, covOp);

  lisU = std::make_shared<Eigen::MatrixXd>(solver.eigenvectors());
  lisW = std::make_shared<Eigen::MatrixXd>(prior->ApplyPrecision(solver.eigenvectors()));

  Eigen::VectorXd deltaVec = solver.eigenvalues().array()/(1.0+solver.eigenvalues().array());
  Eigen::MatrixXd subCov = solver.eigenvectors().transpose() * (*lisW);
  subCov -= deltaVec.asDiagonal();

  lisL = std::make_shared<Eigen::MatrixXd>(subCov.selfadjointView<Eigen::Lower>().llt().matrixL());

  lisToFull = std::make_shared<LIS2Full>(lisU,lisL);
  fullToCS = std::make_shared<CSProjector>(lisU, lisW);

  // Create a new graph for the split likelihood
  std::shared_ptr<WorkGraph> graph = std::make_shared<WorkGraph>();
  graph->AddNode(std::make_shared<SumPiece>(prior->Dimension(), 2), "Parameters");
  graph->AddNode(lisToFull, "Informed Parameters");
  graph->AddNode(fullToCS, "Complementary Parameters");
  graph->AddNode(logLikelihood, "Likelihood");
  graph->AddNode(prior->AsDensity(), "Prior");
  graph->AddNode(std::make_shared<DensityProduct>(2), "Posterior");
  graph->AddEdge("Informed Parameters",0, "Parameters",0);
  graph->AddEdge("Complementary Parameters",0, "Parameters",1);
  graph->AddEdge("Parameters",0,"Prior",0);
  graph->AddEdge("Parameters",0,"Likelihood",0);
  graph->AddEdge("Prior",0,"Posterior",0);
  graph->AddEdge("Likelihood",0,"Posterior",1);

  auto prob = std::make_shared<SamplingProblem>(graph->CreateModPiece("Posterior"));


  lisKernelOpts.put("BlockIndex",0);
  lisKernel = TransitionKernel::Construct(lisKernelOpts,prob);

  csKernelOpts.put("BlockIndex",1);
  csKernel = TransitionKernel::Construct(csKernelOpts,prob);
}

void DILIKernel::UpdateLIS(std::vector<Eigen::VectorXd> const& currState)
{
assert(false);
}
