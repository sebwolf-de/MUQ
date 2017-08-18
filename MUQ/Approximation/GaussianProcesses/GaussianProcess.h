#ifndef GAUSSIANPROCESS_H_
#define GAUSSIANPROCESS_H_

#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"

#include <Eigen/Core>

#include <set>

//#include <nlopt.h>


namespace muq
{

    namespace Approximation
    {


    /** @defgroup GaussianProcesses Gaussian Processes
        @ingroup Approximation
        @brief Tools for defining and working with Gaussian processes.
        @details This group contains classes for defining Gaussian processes and performing Gaussian process regression.

 @code{.cpp}

    const unsigned numFields = 1;   // The number of components in the GP predictions
    const unsigned dim       = 1;   // The number of predictor variables (e.g., location)
    const unsigned numobs    = 100; // The number of observations
    const unsigned numPred   = 10;  // The number of locations where we want to evaluate the GP

    Eigen::MatrixXd trainLocs(dim,       numObs);
    Eigen::MatrixXd trainData(numFields, numObs);

    Eigen::MatrixXd predLocs(dim, numPred);

    // ... Fill in the training locations and data
    // ... Fill in the prediction locations

    // Create the covariance kernel
    const double lengthScale = 0.35;
    const double variance    = 2.0;

    auto kernel = SquaredExpKernel(dim, variance, lengthScael);

    // Create the Mean function
    ConstantMean mean(dim, 1);

    // Create the Gaussian Process
    auto gp = ConstructGP(mean, kernel);

    // Optimize the covariance kernel hyperparameters and store the training data
    gp.Fit(trainLocs, trainData);

    // Make a prediction at new locaitons
    GaussianInformation post = gp.Predict(predLocs);

@endcode

    */

    class GaussianProcess;
    class ObservationInformation;

    struct OptInfo
    {
	GaussianProcess *gp;
    };

    double nlopt_obj(unsigned n, const double *x, double *nlopt_grad, void *opt_info);




    /** @defgroup MeanFunctions Mean Functions
        @ingroup GaussianProcesses
    */
    /** @class MeanFunctionBase
        @ingroup MeanFunctions
    */
    class MeanFunctionBase : public std::enable_shared_from_this<MeanFunctionBase>
    {
    public:
        MeanFunctionBase(unsigned dimIn,
	                 unsigned coDimIn) : inputDim(dimIn), coDim(coDimIn){}

	virtual Eigen::MatrixXd Evaluate(Eigen::MatrixXd const& xs) const = 0;

	virtual std::shared_ptr<MeanFunctionBase> Clone() const = 0;
        virtual std::shared_ptr<MeanFunctionBase> GetPtr()
        {
            return shared_from_this();
        };

	const unsigned inputDim;
	const unsigned coDim;
    };


    /** @class ConstantMean
        @ingroup MeanFunctions
    */
    class ConstantMean : public MeanFunctionBase
    {

    public:
        ConstantMean(unsigned dim, unsigned coDim) : MeanFunctionBase(dim,coDim){};

        virtual std::shared_ptr<MeanFunctionBase> Clone() const override
	{
	    return std::make_shared<ConstantMean>(*this);
	}

	virtual Eigen::MatrixXd Evaluate(Eigen::MatrixXd const& xs) const override
	{
	    return Eigen::MatrixXd::Zero(coDim, xs.cols());
	}
    };

    /** @class LinearTransformMean
        @ingroup MeanFunctions
    */
    template<typename LinearOperator>
    class LinearTransformMean : public MeanFunctionBase
    {

    public:
	template<typename MeanType>
	LinearTransformMean(LinearOperator const& Ain,
			    MeanType const& meanIn) :
	    MeanFunctionBase(meanIn.inputDim, A.rows()),
	    A(Ain),
	    otherMean(meanIn.Clone())
	{
	    assert(A.cols() == otherMean->coDim);
	};

        virtual std::shared_ptr<MeanFunctionBase> Clone() const override
	{
	    return std::make_shared<LinearTransformMean>(*this);
	}

	virtual Eigen::MatrixXd Evaluate(Eigen::MatrixXd const& xs) const override
	{
	    return A * otherMean->Evaluate(xs);
	}

    private:
	LinearOperator A;
	std::shared_ptr<MeanFunctionBase> otherMean;
    };


    template<typename MeanType, typename = typename std::enable_if<std::is_base_of<MeanFunctionBase, MeanType>::value, MeanType>::type>
    LinearTransformMean<MeanType> operator*(Eigen::MatrixXd const& A, MeanType const&K)
    {
	return LinearTransformMean<Eigen::MatrixXd>(A,K);
    }

    /** @class SumMean
        @ingroup MeanFunctions
    */
    class SumMean : public MeanFunctionBase
    {

    public:
	template<typename MeanType1, typename MeanType2>
	SumMean(MeanType1 const& mu1In,
		MeanType2 const& mu2In) :
	    MeanFunctionBase(mu1In.inputDim, mu1In.coDim),
		mu1(mu1In), mu2(mu2In)
	{
	    assert(mu1->inputDim == mu2->inputDim);
	    assert(mu1->coDim == mu2->coDim);
	};

        virtual std::shared_ptr<MeanFunctionBase> Clone() const override
	{
	    return std::make_shared<SumMean>(*this);
	}

	virtual Eigen::MatrixXd Evaluate(Eigen::MatrixXd const& xs) const override
	{
	    return mu1->Evaluate(xs) + mu2->Evaluate(xs);
	}

    private:
	std::shared_ptr<MeanFunctionBase> mu1, mu2;
    };

    template<typename MeanType1, typename MeanType2, typename = typename std::enable_if<std::is_base_of<MeanFunctionBase, MeanType1>::value, MeanType1>::type>
    SumMean operator+(MeanType1 const& mu1, MeanType2 const& mu2)
    {
	return SumMean(mu1, mu2);
    }


    class ObservationInformation
    {
    public:

        ObservationInformation(std::shared_ptr<muq::Utilities::LinearOperator> Hin,
                               Eigen::Ref<const Eigen::VectorXd> const&        locIn,
                               Eigen::Ref<const Eigen::VectorXd> const&        obsIn,
                               Eigen::Ref<const Eigen::MatrixXd> const&        obsCovIn) : H(Hin), loc(locIn), obs(obsIn), obsCov(obsCovIn){};

        // The observation operator
        std::shared_ptr<muq::Utilities::LinearOperator> H;

        // The location of the observation
        Eigen::VectorXd loc;

        // The observed data
        Eigen::VectorXd obs;

        // The covariance of the observational noise
        Eigen::MatrixXd obsCov;
    };



    /** @class GaussianProcess
        @ingroup GaussianProcesses
    */
    class GaussianProcess
    {

    public:

        enum CovarianceType
        {
            DiagonalCov,
            BlockCov,
            FullCov,
            NoCov
        };

        GaussianProcess(MeanFunctionBase& meanIn,
                        KernelBase&       kernelIn) : GaussianProcess(meanIn.Clone(), kernelIn.Clone()){};

        GaussianProcess(std::shared_ptr<MeanFunctionBase> meanIn,
            			std::shared_ptr<KernelBase>       covKernelIn);

        virtual GaussianProcess& Condition(Eigen::Ref<const Eigen::MatrixXd> const& loc,
                                           Eigen::Ref<const Eigen::MatrixXd> const& vals,
                                           double                                   obsVar=0.0);

        virtual GaussianProcess& Condition(std::shared_ptr<ObservationInformation> obs);

        virtual void Optimize();
        
        // Evaluate the mean and covariance
        virtual std::pair<Eigen::MatrixXd, Eigen::MatrixXd> Predict(Eigen::MatrixXd const& newLocs,
                                                                    CovarianceType         covType);

        virtual Eigen::MatrixXd PredictMean(Eigen::MatrixXd const& newPts);

        // Draw a random sample from the GP
        virtual Eigen::MatrixXd Sample(Eigen::MatrixXd const& newPts);


        virtual double LogLikelihood(Eigen::MatrixXd const& xs,
                                     Eigen::MatrixXd const& vals);
        
        // Evaluates the log marginal likelihood needed when fitting hyperparameters
        virtual double MarginalLogLikelihood();
        virtual double MarginalLogLikelihood(Eigen::Ref<Eigen::VectorXd> grad,
                                             bool                        computeGrad=true);

        
      	std::shared_ptr<MeanFunctionBase> Mean(){return mean;};
        std::shared_ptr<KernelBase>       Kernel(){return covKernel;};
        
    protected:

        Eigen::MatrixXd BuildCrossCov(Eigen::MatrixXd const& newLocs);

        void ProcessObservations();


        std::shared_ptr<MeanFunctionBase> mean;
        std::shared_ptr<KernelBase>       covKernel;

        std::vector<std::shared_ptr<ObservationInformation>> observations;


        Eigen::VectorXd trainDiff;
        Eigen::VectorXd sigmaTrainDiff;

        Eigen::LDLT<Eigen::MatrixXd> covSolver;

        int obsDim;
	    const int inputDim;
	    const int coDim;

        // Have new observations been added since the covariance was inverted?
        bool hasNewObs;

	    const double pi = 4.0 * atan(1.0); //boost::math::constants::pi<double>();

    };// class GaussianProcess

    /** @ingroup GaussianProcesses
     */
    template<typename MeanType, typename KernelType>
    GaussianProcess ConstructGP(MeanType const& mean,
				KernelType const& kernel)
    {
	return GaussianProcess(mean.Clone(),kernel.Clone());
    }


    } // namespace Approximation
} // namespace muq

#endif // #ifndef GAUSSIANPROCESS_H_
