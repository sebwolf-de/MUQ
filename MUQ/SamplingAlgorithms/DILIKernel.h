#ifndef DILIKERNEL_H_
#define DILIKERNEL_H_

#include "MUQ/SamplingAlgorithms/TransitionKernel.h"

#include "MUQ/SamplingAlgorithms/MCMCProposal.h"
#include "MUQ/Modeling/Distributions/GaussianBase.h"
#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

namespace muq {
  namespace SamplingAlgorithms {

    class CSProjector : public muq::Modeling::LinearOperator
    {
    public:

      CSProjector(std::shared_ptr<Eigen::MatrixXd> const& Uin,
                  std::shared_ptr<Eigen::MatrixXd> const& Win) : LinearOperator(Uin->rows(), Win->rows()),
                                                                 U(Uin), W(Win){};

      virtual ~CSProjector() = default;

      /** Apply the linear operator to a vector */
      virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

      /** Apply the transpose of the linear operator to a vector. */
      virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

    private:
      std::shared_ptr<Eigen::MatrixXd> U, W;
    };


    class LIS2Full : public muq::Modeling::LinearOperator
    {
    public:
      LIS2Full(std::shared_ptr<Eigen::MatrixXd> const& Uin,
               std::shared_ptr<Eigen::MatrixXd> const& Lin) : LinearOperator(Uin->rows(), Uin->cols()),
                                                              U(Uin), L(Lin){};

      virtual ~LIS2Full() = default;

      /** Apply the linear operator to a vector */
      virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

      /** Apply the transpose of the linear operator to a vector. */
      virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

    private:
      std::shared_ptr<Eigen::MatrixXd> U, L;
    };

    /**
      @ingroup MCMCKernels
      @class DILIKernel
      @brief An implementation of the Dimension Independent Sikelihood Informed subspace (DILI) MCMC sampler
      @details <B>Configuration Parameters:</B>

      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "LIS Kernel"  | String | - | A string pointing to a block of kernel/proposal options for the Likelihood informed subspace. |
      "CS Kernel"   | String | - | A string pointing to a block of kernel/proposal options for the Complementary space.  Typically this will be a Crank-Nicolson proposal. |
     */
    class DILIKernel : public TransitionKernel {
    public:

      DILIKernel(boost::property_tree::ptree       const& pt,
                 std::shared_ptr<AbstractSamplingProblem> problem);


      DILIKernel(boost::property_tree::ptree                  const& pt,
                 std::shared_ptr<AbstractSamplingProblem>            problem,
                 std::shared_ptr<muq::Modeling::ModPiece>     const& likelihood,
                 std::shared_ptr<muq::Modeling::GaussianBase> const& prior);

      virtual ~DILIKernel() = default;

      virtual inline std::shared_ptr<TransitionKernel> LISKernel() {return lisKernel;};
      virtual inline std::shared_ptr<TransitionKernel> CSKernel() {return csKernel;};

      virtual void PostStep(unsigned int const t,
                            std::vector<std::shared_ptr<SamplingState>> const& state) override;

      virtual std::vector<std::shared_ptr<SamplingState>> Step(unsigned int const t,
                                                               std::shared_ptr<SamplingState> prevState) override;

      virtual void PrintStatus(std::string prefix) const override;

    protected:

      /** From a ModGraphPiece defining the posterior log density, this function
          extracts a ModGraphPiece defining the likelihood function.
      */
      static std::shared_ptr<muq::Modeling::ModPiece> ExtractLikelihood(std::shared_ptr<muq::Modeling::ModPiece> const& tgtGraphPiece,
                                                                        std::string                              const& nodeName="Likelihood");


      /** From a ModGraphPiece defining the posterior log density, this function
          extracts a GaussianBase instance defining the prior.
      */
      static std::shared_ptr<muq::Modeling::GaussianBase> ExtractPrior(std::shared_ptr<muq::Modeling::ModPiece> const& tgtGraphPiece,
                                                                       std::string                              const& nodeName="Prior");

      /** Create the likelihood informed subspace for the first time. */
      void CreateLIS(std::vector<Eigen::VectorXd> const& currState,
                     unsigned int                        blockId,
                     boost::property_tree::ptree  const& solverOpts);

      /** Update the likelihood informed subspace using Hessian information at a
          new point.
      */
      void UpdateLIS(std::vector<Eigen::VectorXd> const& currState);

      boost::property_tree::ptree lisKernelOpts;
      boost::property_tree::ptree csKernelOpts;

      std::shared_ptr<muq::Modeling::ModPiece> logLikelihood;
      std::shared_ptr<muq::Modeling::GaussianBase> prior;

      std::shared_ptr<muq::Modeling::ModPiece> forwardModel;
      std::shared_ptr<muq::Modeling::GaussianBase> noiseModel;

      // A matrix containing the eigenvectors of the generalized eigenvalue problem Hv = lam*\Gamma^{-1}v
      std::shared_ptr<Eigen::MatrixXd> lisU;

      // W = \Gamma^{-1} v
      std::shared_ptr<Eigen::MatrixXd> lisW;

      // L is the Cholesky factor of the approximate posterior covariance projected onto the subspace
      std::shared_ptr<Eigen::MatrixXd> lisL;

      // Defines the operation from the LIS to the full space
      std::shared_ptr<muq::Modeling::LinearOperator> lisToFull;

      // Defines the projection onto the complementary space
      std::shared_ptr<muq::Modeling::LinearOperator> fullToCS;

      // Transition kernel for updating the LIS
      std::shared_ptr<TransitionKernel> lisKernel;

      // Transition kernel for updating the CS
      std::shared_ptr<TransitionKernel> csKernel;

      // The type of Hessian approximation to use.  "Exact" or "GaussNewton"
      std::string hessType;

      // Options for the Generalized eigenvalue solver
      boost::property_tree::ptree lobpcgOpts;

    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
