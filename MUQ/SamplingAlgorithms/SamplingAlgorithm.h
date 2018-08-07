#ifndef SAMPLINGALGORITHM_H_
#define SAMPLINGALGORITHM_H_

#include "MUQ/Modeling/WorkPiece.h"

#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SampleCollection.h"

namespace muq {
  namespace SamplingAlgorithms {

    class SamplingAlgorithm {//} : public muq::Modeling::WorkPiece {
    public:

      SamplingAlgorithm(std::shared_ptr<SampleCollection> samplesIn,
                        std::shared_ptr<SampleCollection> QOIsIn)
       : samples(samplesIn),
         QOIs(QOIsIn)
       {}

      SamplingAlgorithm(std::shared_ptr<SampleCollection> samplesIn)
       : SamplingAlgorithm (samplesIn, std::make_shared<SampleCollection>())
      {}

      virtual ~SamplingAlgorithm() = default;

      std::shared_ptr<SampleCollection> GetSamples() const{return samples;};
      std::shared_ptr<SampleCollection> GetQOIs() const{return QOIs;};

      virtual std::shared_ptr<SampleCollection> Run(){return RunImpl();};


      virtual std::shared_ptr<SampleCollection> RunImpl() = 0;

    protected:

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

      std::shared_ptr<SampleCollection> QOIs;

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
