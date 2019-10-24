#ifndef PYDISTRIBUTION_H_
#define PYDISTRIBUTION_H_

#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/Modeling/Distributions/GaussianBase.h"


namespace muq {
  namespace Modeling {
    class PyDistribution : public Distribution {
    public:
      PyDistribution(unsigned int varSizeIn, Eigen::VectorXi const& hyperSizesIn = Eigen::VectorXi());

      static std::vector<Eigen::VectorXd> ToStdVec(ref_vector<Eigen::VectorXd> const& input);

    protected:

      virtual Eigen::VectorXd SampleImpl(std::vector<Eigen::VectorXd> const& inputs) = 0;
      virtual Eigen::VectorXd SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) override;


      virtual double LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual double LogDensityImpl(std::vector<Eigen::VectorXd> const& inputs) = 0;
    };


    class PyGaussianBase : public GaussianBase {
    public:

      using GaussianBase::GaussianBase;

      virtual ~PyGaussianBase() = default;

      virtual Eigen::MatrixXd ApplyCovariance(Eigen::Ref<const Eigen::MatrixXd> const& x) const override;
      virtual Eigen::MatrixXd ApplyCovariance(Eigen::MatrixXd const& x) const = 0;

      virtual Eigen::MatrixXd ApplyPrecision(Eigen::Ref<const Eigen::MatrixXd> const& x) const override;
      virtual Eigen::MatrixXd ApplyPrecision(Eigen::MatrixXd const& x) const = 0;

      virtual Eigen::MatrixXd ApplyCovSqrt(Eigen::Ref<const Eigen::MatrixXd> const& x) const override;
      virtual Eigen::MatrixXd ApplyCovSqrt(Eigen::MatrixXd const& x) const = 0;

      virtual Eigen::MatrixXd ApplyPrecSqrt(Eigen::Ref<const Eigen::MatrixXd> const& x) const override;
      virtual Eigen::MatrixXd ApplyPrecSqrt(Eigen::MatrixXd const& x) const = 0;

      virtual void ResetHyperparameters(ref_vector<Eigen::VectorXd> const& params) override;
      virtual void ResetHyperparameters(std::vector<Eigen::VectorXd> const& params){};
    };

  } // namespace Modeling
} // namespace muq

#endif
