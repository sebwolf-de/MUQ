#include "MUQ/Approximation/GaussianProcesses/ObservationInformation.h"


using namespace muq::Utilities;
using namespace muq::Approximation;


void ObservationInformation::FillCrossCov(std::shared_ptr<ObservationInformation> otherObs,
                                          std::shared_ptr<KernelBase>              kernel,
                                          Eigen::Ref<Eigen::MatrixXd>             covBlock)
{
    // check to see if the other observation is a derivative observation
    std::shared_ptr<DerivativeObservation> derivObs = std::dynamic_pointer_cast<DerivativeObservation>(otherObs);
    if(derivObs){
        kernel->FillDerivCovariance(loc, derivObs->loc, derivObs->derivCoords, covBlock);
    }else{ 
        kernel->FillCovariance(loc, otherObs->loc, covBlock);
    }
};

void DerivativeObservation::FillCrossCov(std::shared_ptr<ObservationInformation> otherObs,
                                         std::shared_ptr<KernelBase>              kernel,
                                         Eigen::Ref<Eigen::MatrixXd>             covBlock)
{
    // check to see if the other observation is a derivative observation
    std::shared_ptr<DerivativeObservation> derivObs = std::dynamic_pointer_cast<DerivativeObservation>(otherObs);
    if(derivObs){
        // concatenate the derivatives together
        std::vector<unsigned> newDerivCoords = derivCoords;
        newDerivCoords.insert(newDerivCoords.end(), derivObs->derivCoords.begin(), derivObs->derivCoords.end());
        kernel->FillDerivCovariance(loc, otherObs->loc, newDerivCoords, covBlock);
    }else{
        kernel->FillDerivCovariance(loc, otherObs->loc, derivCoords, covBlock);
    }
};

