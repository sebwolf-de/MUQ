#ifndef SAMPLINGALGORITHM_H_
#define SAMPLINGALGORITHM_H_

#include "MUQ/Modeling/WorkPiece.h"

#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SampleCollection.h"

namespace muq {
  namespace SamplingAlgorithms {

    class SamplingAlgorithm {//} : public muq::Modeling::WorkPiece {
    public:

      SamplingAlgorithm(std::shared_ptr<SampleCollection> samplesIn) : samples(samplesIn){};

      virtual ~SamplingAlgorithm() = default;

      std::shared_ptr<SampleCollection> GetSamples() const{return samples;};

      virtual std::shared_ptr<SampleCollection> Run(){return Run(std::vector<Eigen::VectorXd>());};
      virtual std::shared_ptr<SampleCollection> Run(Eigen::VectorXd const& x0){return Run(std::vector<Eigen::VectorXd>(1,x0));};
      virtual std::shared_ptr<SampleCollection> Run(std::vector<Eigen::VectorXd> const& x0){ return RunImpl(x0);};

    protected:

      virtual std::shared_ptr<SampleCollection> RunImpl(std::vector<Eigen::VectorXd> const& x0) = 0;

      /**
	 Inputs:
	 <ol>
	 <li> Parameters for the algorithm
	 <li> The muq::SamplingAlgorithms::SamplingProblem that evaluates/samples the target distribution
	 </ol>
	 @param[in] inputs Inputs to the algorithm
       */
      //virtual void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      std::shared_ptr<SampleCollection> samples;

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
