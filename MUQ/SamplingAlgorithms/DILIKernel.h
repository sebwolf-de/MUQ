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
    @class AverageHessian
    @details During the construction of the LIS, we need the average Hessian.
             This class can be used to compute the action of the average Hessian
             given a generalized eigenvalue decomposition of the old Hessian and
             a linear operator for the new Hessian, i.e.,
             $$
             \bar{H} = \frac{m}{m+1}H_{old} + \frac{1}{m+1}H_{new}
             $$
    */
    class AverageHessian : public muq::Modeling::LinearOperator
    {
    public:
      AverageHessian(unsigned int                            numOldSamps,
                     std::shared_ptr<Eigen::MatrixXd> const& oldUIn,
                     std::shared_ptr<Eigen::MatrixXd> const& oldWIn,
                     std::shared_ptr<Eigen::VectorXd> const& oldValsIn,
                     std::shared_ptr<muq::Modeling::LinearOperator> const& newHess);

      /** Apply the linear operator to a vector */
      virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

      /** Apply the transpose of the linear operator to a vector. */
      virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) override;

    private:
      const double numSamps;
      std::shared_ptr<Eigen::MatrixXd> oldU, oldW;
      std::shared_ptr<Eigen::VectorXd> oldEigVals;
      std::shared_ptr<muq::Modeling::LinearOperator> newHess;

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
                 std::shared_ptr<muq::Modeling::GaussianBase> const& prior,
                 std::shared_ptr<muq::Modeling::ModPiece>     const& noiseModel,
                 std::shared_ptr<muq::Modeling::ModPiece>     const& forwardModelIn);

      DILIKernel(boost::property_tree::ptree                  const& pt,
                 std::shared_ptr<AbstractSamplingProblem>            problem,
                 std::shared_ptr<muq::Modeling::GaussianBase> const& prior,
                 std::shared_ptr<muq::Modeling::ModPiece>     const& likelihood);

      virtual ~DILIKernel() = default;

      virtual inline std::shared_ptr<TransitionKernel> LISKernel() {return lisKernel;};
      virtual inline std::shared_ptr<TransitionKernel> CSKernel() {return csKernel;};

      virtual void PostStep(unsigned int const t,
                            std::vector<std::shared_ptr<SamplingState>> const& state) override;

      virtual std::vector<std::shared_ptr<SamplingState>> Step(unsigned int const t,
                                                               std::shared_ptr<SamplingState> prevState) override;

      virtual void PrintStatus(std::string prefix) const override;

      /** From a ModGraphPiece defining the posterior log density, this function
          extracts a ModGraphPiece defining the likelihood function.
      */
      static std::shared_ptr<muq::Modeling::ModPiece> ExtractLikelihood(std::shared_ptr<AbstractSamplingProblem> const& problem,
                                                                        std::string                              const& nodeName);


      /** From a ModGraphPiece defining the posterior log density, this function
          extracts a GaussianBase instance defining the prior.
      */
      static std::shared_ptr<muq::Modeling::GaussianBase> ExtractPrior(std::shared_ptr<AbstractSamplingProblem> const& problem,
                                                                       std::string                              const& nodeName);

      static std::shared_ptr<muq::Modeling::ModPiece> ExtractNoiseModel(std::shared_ptr<muq::Modeling::ModPiece> const& likelihood);

      // Extract the forward model from the likelihood.  This function assumes the graph has only one node between the output of the likelihood function and the output of the forward model
      static std::shared_ptr<muq::Modeling::ModPiece> ExtractForwardModel(std::shared_ptr<muq::Modeling::ModPiece> const& likelihoodIn);

      static std::shared_ptr<muq::Modeling::ModPiece> CreateLikelihood(std::shared_ptr<muq::Modeling::ModPiece> const& forwardModel,
                                                                       std::shared_ptr<muq::Modeling::ModPiece> const& noiseDensity);


    protected:

      /** Create the likelihood informed subspace for the first time. */
      void CreateLIS(std::vector<Eigen::VectorXd> const& currState);

      /** Update the likelihood informed subspace using Hessian information at a
          new point.
          @param[in] numSamps The number of samples that have been used to estimate the current LIS.
          @param[in] currState The current state of the chain.
      */
      void UpdateLIS(unsigned int                        numSamps,
                     std::vector<Eigen::VectorXd> const& currState);

      // Used to reconstruct the LIS and CS kernels when lisU, lisW, or lisL changes.
      void UpdateKernels();

      boost::property_tree::ptree lisKernelOpts;
      boost::property_tree::ptree csKernelOpts;

      std::shared_ptr<muq::Modeling::ModPiece> logLikelihood;
      std::shared_ptr<muq::Modeling::GaussianBase> prior;

      std::shared_ptr<muq::Modeling::ModPiece> forwardModel;
      std::shared_ptr<muq::Modeling::ModPiece> noiseDensity;

      // A matrix containing the eigenvectors of the generalized eigenvalue problem Hv = lam*\Gamma^{-1}v
      std::shared_ptr<Eigen::MatrixXd> lisU;

      // A vector of the LIS eigenvalues
      std::shared_ptr<Eigen::VectorXd> lisEigVals;

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

      const unsigned int updateInterval;
      unsigned int numLisUpdates;
    };
  } // namespace SamplingAlgorithms
} // namespace muq

#endif
