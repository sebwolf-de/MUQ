
class MCSampleProposal : public MCMCProposal {
public:
  MCSampleProposal(boost::property_tree::ptree       const& pt,
                   std::shared_ptr<AbstractSamplingProblem> prob,
                   std::shared_ptr<Distribution> dist
                  )
   : MCMCProposal(pt, prob),
     dist(dist)
  {}

  std::shared_ptr<SamplingState> Sample(std::shared_ptr<SamplingState> const& currentState) override {
    return std::make_shared<SamplingState>(dist->Sample());
  }

  double LogDensity(std::shared_ptr<SamplingState> const& currState,
                    std::shared_ptr<SamplingState> const& propState) override {
    return 0.0;
  }


private:
  std::shared_ptr<Distribution> dist;
};
