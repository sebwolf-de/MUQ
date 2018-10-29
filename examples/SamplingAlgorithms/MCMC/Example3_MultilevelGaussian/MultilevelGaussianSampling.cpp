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






class MySamplingProblem : public AbstractSamplingProblem {
public:
  MySamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> targetIn)
   : AbstractSamplingProblem(Eigen::VectorXi::Constant(1,2), Eigen::VectorXi::Constant(1,2)),
     target(targetIn){}

  virtual ~MySamplingProblem() = default;


  virtual double LogDensity(std::shared_ptr<SamplingState> state) override {
    lastState = state;
    return target->Evaluate(state->state).at(0)(0);
  };

  virtual Eigen::VectorXd GradLogDensity(std::shared_ptr<SamplingState> state,
                                         unsigned                       blockWrt) override {
    assert (false);
  };

  virtual std::shared_ptr<SamplingState> QOI() override {
    assert (lastState != nullptr);
    return std::make_shared<SamplingState>(lastState->state, 1.0);
  }

private:
  std::shared_ptr<SamplingState> lastState = nullptr;

  std::shared_ptr<muq::Modeling::ModPiece> target;

};


class MyInterpolation : public MIInterpolation {
public:
  std::shared_ptr<SamplingState> interpolate (std::shared_ptr<SamplingState> coarseProposal, std::shared_ptr<SamplingState> fineProposal) {
    return std::make_shared<SamplingState>(coarseProposal->state);
  }
};

class MyMIComponentFactory : public MIComponentFactory {
public:
  virtual std::shared_ptr<MCMCProposal> proposal (std::shared_ptr<MultiIndex> index, std::shared_ptr<AbstractSamplingProblem> samplingProblem) override {
    pt::ptree pt;
    pt.put("BlockIndex",0);

    Eigen::VectorXd mu(2);
    mu << 1.0, 2.0;
    Eigen::MatrixXd cov(2,2);
    cov << 0.7, 0.6,
    0.6, 1.0;
    cov *= 20.0;

    auto prior = std::make_shared<Gaussian>(mu, cov);

    return std::make_shared<CrankNicolsonProposal>(pt, samplingProblem, prior);
  }

  virtual std::shared_ptr<MultiIndex> finestIndex() override {
    auto index = std::make_shared<MultiIndex>(1);
    index->SetValue(0, 3);
    return index;
  }

  virtual std::shared_ptr<MCMCProposal> coarseProposal (std::shared_ptr<MultiIndex> index,
                                                        std::shared_ptr<AbstractSamplingProblem> coarseProblem,
                                                           std::shared_ptr<SingleChainMCMC> coarseChain) override {
    pt::ptree ptProposal;
    ptProposal.put("BlockIndex",0);
    int subsampling = 5;
    ptProposal.put("subsampling", subsampling);
    return std::make_shared<SubsamplingMIProposal>(ptProposal, coarseProblem, coarseChain);
  }

  virtual std::shared_ptr<AbstractSamplingProblem> samplingProblem (std::shared_ptr<MultiIndex> index) override {
    Eigen::VectorXd mu(2);
    mu << 1.0, 2.0;
    Eigen::MatrixXd cov(2,2);
    cov << 0.7, 0.6,
           0.6, 1.0;

    if (index->GetValue(0) == 0) {
      mu *= 0.8;
      cov *= 2.0;
    } else if (index->GetValue(0) == 1) {
      mu *= 0.9;
      cov *= 1.5;
    } else if (index->GetValue(0) == 2) {
      mu *= 0.99;
      cov *= 1.1;
    } else if (index->GetValue(0) == 3) {
      mu *= 1.0;
      cov *= 1.0;
    } else {
      std::cerr << "Sampling problem not defined!" << std::endl;
      assert (false);
    }

    auto coarseTargetDensity = std::make_shared<Gaussian>(mu, cov)->AsDensity();
    return std::make_shared<MySamplingProblem>(coarseTargetDensity);
  }

  virtual std::shared_ptr<MIInterpolation> interpolation (std::shared_ptr<MultiIndex> index) override {
    return std::make_shared<MyInterpolation>();
  }

  virtual Eigen::VectorXd startingPoint (std::shared_ptr<MultiIndex> index) override {
    Eigen::VectorXd mu(2);
    mu << 1.0, 2.0;
    return mu;
  }

};

int main(){

  auto componentFactory = std::make_shared<MyMIComponentFactory>();

  pt::ptree pt;

  pt.put("NumSamples", 1e4); // number of samples for single level
  pt.put("NumInitialSamples", 1e3); // number of initial samples for greedy MLMCMC
  pt.put("GreedyTargetVariance", 0.05); // estimator variance to be achieved by greedy algorithm

  std::cout << std::endl << "*************** greedy multillevel chain" << std::endl << std::endl;

  GreedyMLMCMC greedymlmcmc (pt, componentFactory);
  greedymlmcmc.run();
  std::cout << "mean QOI: " << greedymlmcmc.meanQOI().transpose() << std::endl;
  greedymlmcmc.draw(false);

  std::cout << std::endl << "*************** single chain reference" << std::endl << std::endl;

  SLMCMC slmcmc (pt, componentFactory);
  slmcmc.run();
  std::cout << "mean QOI: " << slmcmc.meanQOI().transpose() << std::endl;

  return 0;
}
