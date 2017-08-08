#ifndef KARHUNENLOEVEEXPANSION_H
#define KARHUNENLOEVEEXPANSION_H

#include <Eigen/Core>

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"

#include <boost/property_tree/ptree.hpp>

namespace muq
{
namespace Approximation
{

    class KarhunenLoeveExpansion
    {

    public:

        KarhunenLoeveExpansion(std::shared_ptr<KernelBase> kernelIn,
                               Eigen::MatrixXd      const& seedPtsIn,
                               Eigen::VectorXd      const& seedWtsIn,
                               boost::property_tree::ptree options = boost::property_tree::ptree());

        /** Evaluates the KL modes at one or more locations.  Each column of the pts matrix contains a point where we want to evaluate the modes.  Each column of the output contains a mode.  Each row of the output corresponds to an input point.
        */
        Eigen::MatrixXd GetModes(Eigen::Ref<const Eigen::MatrixXd> const& pts);


        Eigen::MatrixXd Evaluate(Eigen::Ref<const Eigen::MatrixXd> const& pt,
                                 Eigen::Ref<const Eigen::VectorXd> const& coeffs);

    private:


        // Points used to discretize the KL modes.
        Eigen::MatrixXd seedPts;
        Eigen::VectorXd seedWts;

        // The covariance kernel used to construct this expansion
        std::shared_ptr<KernelBase> covKernel;

        // Values of the KL modes at the seed points.  Each column corresponds to a basis function and each row to a pt function
        Eigen::MatrixXd modeVecs;
        Eigen::VectorXd modeEigs;


    }; // class KarhuneLoeveExpansion



}// namespace Approximation
}// namespace muq


#endif
