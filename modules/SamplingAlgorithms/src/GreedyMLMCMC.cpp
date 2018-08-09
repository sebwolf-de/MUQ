#include "MUQ/SamplingAlgorithms/GreedyMLMCMC.h"

namespace muq {
  namespace SamplingAlgorithms {

    GreedyMLMCMC::GreedyMLMCMC (pt::ptree pt, std::shared_ptr<MIComponentFactory> componentFactory)
    : componentFactory(componentFactory),
      numInitialSamples(pt.get("NumInitialSamples",1000))
    {

      const int levels = componentFactory->finestIndex()->GetValue(0);
      const auto rootIndex = std::make_shared<MultiIndex>(1);


      const double e = pt.get("GreedyTargetVariance",0.1);
      const double beta = pt.get("GreedyResamplingFactor",0.5);


      for (int level = 0; level <= levels; level++) {

        auto boxHighestIndex = std::make_shared<MultiIndex>(1,level);
        auto box = std::make_shared<MIMCMCBox>(componentFactory, boxHighestIndex);
        boxes.push_back(box);

        for (int samp = 0; samp < numInitialSamples; samp++) {
          box->Sample();
        }
      }


      while(true) {
        double var_mle = 0.0;
        for (int i = 0; i <= levels; i++) {
          auto qois = boxes[i]->finestChain()->GetQOIs();
          var_mle += qois->Variance().cwiseQuotient(qois->ESS()).maxCoeff();
        }
        if (var_mle <= std::pow(e,2))
          break;

        // Find level with largest payoff ratio
        int l = -1;
        double payoff_l = -1;
        for (int i = 0; i <= levels; i++) {
          auto qois = boxes[i]->finestChain()->GetQOIs();

          double my_payoff = qois->Variance().maxCoeff() / boxes[i]->finestChain()->TotalTime();

          if (my_payoff > payoff_l) {
            l = i;
            payoff_l = my_payoff;
          }
        }

        // Beta percent new samples on largest payoff level
        double weight_sum = 0.0;
        auto finestChain = boxes[l]->finestChain();
        for (int s = 0; s < finestChain->GetSamples()->samples.size(); s++) {
          std::shared_ptr<SamplingState> sample = finestChain->GetSamples()->samples[s];
          weight_sum += sample->weight;
        }
        int n_samples = std::ceil(weight_sum);
        int n_new_samples = std::ceil(n_samples * beta);

        std::cout << "var_mle " << var_mle << "\t" << n_new_samples << " new samples on level " << l << std::endl;
        for (int i = 0; i < n_new_samples; i++)
          boxes[l]->Sample();
      }

      // Compute full QOI estimate
      Eigen::VectorXd MImean(componentFactory->samplingProblem(componentFactory->finestIndex())->blockSizesQOI.sum());
      MImean.setZero();

      for (auto box : boxes) {
        Eigen::VectorXd sampMean = box->Mean();

        std::cout << "Sample Mean = \n" << sampMean.transpose() << std::endl;
        MImean += sampMean;
      }

      std::cout << "ML Mean = \n" << MImean.transpose() << std::endl;
    }

    void GreedyMLMCMC::draw() {
      std::ofstream graphfile;
      graphfile.open ("graph");
      graphfile << "digraph {" << std::endl;
      graphfile << "nodesep=1.2;" << std::endl;
      graphfile << "splines=false;" << std::endl;
      for (auto box : boxes) {
        box->draw(graphfile);
      }
      graphfile << "}" << std::endl;
      graphfile.close();
    }

  }
}
