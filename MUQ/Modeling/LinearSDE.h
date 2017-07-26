#ifndef LINEARSDE_H
#define LINEARSDE_H


#include "MUQ/Utilities/LinearAlgebra/LinearOperator.h"
#include "MUQ/Utilities/LinearAlgebra/EigenLinearOperator.h"

#include <boost/property_tree/ptree.hpp>
#include <random>

namespace muq
{
namespace Modeling
{

    /** @brief Defines a linear time invariant stochastic differential equation.
        @details This class defines a LTI SDE of the form
\f[
\frac{\partial f(t)}{\partial t} = F f(t) + L w(t),
\f]
where \f$f(t)\f$ is the solution in \f$\mathbb{R}^M\f$, \f$F\f$ is an \f$M\timesM\f$ matrix, \f$L\f$ is an \f$M\times N\f$ matrix, and \f$w(t)\f$ is a white noise process with an \f$N\timesN\f$ covariance matrix \f$Q\f$.
    */
    class LinearSDE
    {

    public:

        template<typename Derived1, typename Derived2>
        LinearSDE(Eigen::Matrix<Derived1,Eigen::Dynamic, Eigen::Dynamic>   const& Fin,
                  Eigen::Matrix<Derived2,Eigen::Dynamic, Eigen::Dynamic>   const& Lin,
                  Eigen::MatrixXd                                          const& Qin,
                  boost::property_tree::ptree options) : LinearSDE(muq::Utilities::LinearOperator::Create(Fin),
                                                                   muq::Utilities::LinearOperator::Create(Lin),
                                                                   Qin,
                                                                   options)
        {};
        
        
        LinearSDE(std::shared_ptr<muq::Utilities::LinearOperator>    Fin,
                  std::shared_ptr<muq::Utilities::LinearOperator>    Lin,
                  Eigen::MatrixXd                             const& Qin,
                  boost::property_tree::ptree                        options);


        /** Given \f$f(t)\f$, the state of the system at time \f$t\f$, return a random realization of the state at time \f$t+\delta t\f$.
         */
        Eigen::VectorXd EvolveState(Eigen::VectorXd const& f0,
                                    double                 T) const;

        /** Given the mean and covariance of the solution at time \f$t\f$, compute the mean and covariance of the solution at time \f$t+T\f$.
         */
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> EvolveDistribution(Eigen::VectorXd const& muIn,
                                                                       Eigen::MatrixXd const& gammaIn,
                                                                       double                 T) const; 

        /** Evolve the mean and covariance of the system using a std::pair to hold the distribution.
         */
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> EvolveDistribution(std::pair<Eigen::VectorXd,Eigen::MatrixXd> const& muCov,
                                                                       double                                            T) const{
            return EvolveDistribution(muCov.first, muCov.second, dt);
        }; 


        /// The dimension of the state variable \f$f(t)\f$.
        const int stateDim;


        std::shared_ptr<muq::Utilities::LinearOperator> GetF() const{return F;};
        std::shared_ptr<muq::Utilities::LinearOperator> GetL() const{return L;};
        Eigen::MatrixXd const& GetQ() const{return Q;};
        
        
    protected:

        void ExtractOptions(boost::property_tree::ptree options);
        
        std::shared_ptr<muq::Utilities::LinearOperator> F;
        std::shared_ptr<muq::Utilities::LinearOperator> L;

        Eigen::MatrixXd Q;
        Eigen::MatrixXd sqrtQ;

        double dt; // time step used in SDE integration

    };


}// namespace Modeling
}// namespace muq




#endif
