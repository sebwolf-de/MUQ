#ifndef SAMPLINGALGORITHM_H_
#define SAMPLINGALGORITHM_H_

#include "MUQ/Modeling/WorkPiece.h"

#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SampleCollection.h"

namespace muq {
  namespace SamplingAlgorithms {

    class SamplingAlgorithm {//} : public muq::Modeling::WorkPiece {
    public:

      SamplingAlgorithm(){};

      virtual ~SamplingAlgorithm() = default;

      SampleCollection const& GetSamples() const;

      virtual SampleCollection const& Run(){return Run(std::vector<Eigen::VectorXd>());};
      virtual SampleCollection const& Run(Eigen::VectorXd const& x0){return Run(std::vector<Eigen::VectorXd>(1,x0));};
      virtual SampleCollection const& Run(std::vector<Eigen::VectorXd> const& x0){return RunImpl(x0);};

      virtual SampleCollection const& RunImpl(std::vector<Eigen::VectorXd> const& x0) = 0;

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

      SampleCollection samples;

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
