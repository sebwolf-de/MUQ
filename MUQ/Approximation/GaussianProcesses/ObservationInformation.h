#ifndef OBSERVATIONINFORMATION_H
#define OBSERVATIONINFORMATION_H

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"
#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

#include <Eigen/Core>

#include <vector>

namespace muq{
namespace Approximation{


    class ObservationInformation
    {
    public:

        ObservationInformation(std::shared_ptr<muq::Modeling::LinearOperator> Hin,
                               Eigen::Ref<const Eigen::VectorXd> const&        locIn,
                               Eigen::Ref<const Eigen::VectorXd> const&        obsIn,
                               Eigen::Ref<const Eigen::MatrixXd> const&        obsCovIn) : H(Hin), loc(locIn), obs(obsIn), obsCov(obsCovIn){};

        virtual ~ObservationInformation() = default;

        virtual void FillCrossCov(std::shared_ptr<ObservationInformation> otherObs,
                                  std::shared_ptr<KernelBase>             kernel,
                                  Eigen::Ref<Eigen::MatrixXd>             covBlock);

        // The observation operator
        std::shared_ptr<muq::Modeling::LinearOperator> H;

        // The location of the observation
        Eigen::VectorXd loc;

        // The observed data
        Eigen::VectorXd obs;

        // The covariance of the observational noise
        Eigen::MatrixXd obsCov;

    };


    class DerivativeObservation : public ObservationInformation
    {
    public:

        DerivativeObservation(std::shared_ptr<muq::Modeling::LinearOperator> Hin,
                              Eigen::Ref<const Eigen::VectorXd> const&        locIn,
                              Eigen::Ref<const Eigen::VectorXd> const&        obsIn,
                              Eigen::Ref<const Eigen::MatrixXd> const&        obsCovIn,
                              std::vector<int>                                derivCoordsIn) : ObservationInformation(Hin, locIn, obsIn, obsCovIn), derivCoords(derivCoordsIn){};

        virtual ~DerivativeObservation() = default;

        virtual void FillCrossCov(std::shared_ptr<ObservationInformation> otherObs,
                                  std::shared_ptr<KernelBase>             kernel,
                                  Eigen::Ref<Eigen::MatrixXd>             covBlock) override;

        // coordinate directions of the derivatives.  derivCoords.size() corresponds to the derivative order.  For example, a second derivative will have derivCoords.size()==2.
        std::vector<int> derivCoords;

    };


} // namespace muq
} // namespace Approximation





#endif // #ifndef OBSERVATIONINFORMATION_H
