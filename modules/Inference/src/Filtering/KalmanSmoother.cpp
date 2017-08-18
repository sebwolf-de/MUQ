#include "MUQ/Inference/Filtering/KalmanSmoother.h"

#include <Eigen/Dense>

using namespace muq::Utilities;
using namespace muq::Inference;




std::pair<Eigen::VectorXd, Eigen::MatrixXd> KalmanSmoother::Analyze(std::pair<Eigen::VectorXd, Eigen::MatrixXd> const& currDist_t,
                                                                    std::pair<Eigen::VectorXd, Eigen::MatrixXd> const& nextDist_t,
                                                                    std::pair<Eigen::VectorXd, Eigen::MatrixXd> const& nextDist_n,
                                                                    std::shared_ptr<muq::Utilities::LinearOperator>    F)
{

    std::pair<Eigen::VectorXd, Eigen::MatrixXd> output;

    Eigen::MatrixXd C = nextDist_t.second.llt().solve( F->Apply(currDist_t.second) ).transpose();

    output.first = currDist_t.first + C*(nextDist_n.first - nextDist_t.first);
    output.second = currDist_t.second + C*(nextDist_n.second - nextDist_t.second).selfadjointView<Eigen::Lower>()*C.transpose();

    return output;
}
               
