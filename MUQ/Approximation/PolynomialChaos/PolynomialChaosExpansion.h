#ifndef POLYNOMIALCHAOSEXPANSION_H_
#define POLYNOMIALCHAOSEXPANSION_H_

#include "MUQ/Approximation/Polynomials/BasisExpansion.h"
#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

#include <vector>
#include <set>

namespace muq {
namespace Approximation {

  class PCEFactory;

  /**
   @class PolynomialChaosExpansion
   @ingroup PolynomialChaos
   @brief A class for representing and using expansions of orthogonal multivariate polynomials
   @details
   * A particular polynomial chaos expansion for a function from R^n->R^m. This class uses
   * some the MultiIndexSet class to define PCE terms that are in the expansion.
   * For each PCE term and output, there is a coefficient. The PCE is built from 1D polynomials
   * specified so each dimension may have independent polynomial
   * choices.
   *
  @see BasisExpansion, PCEFactory
   */
  class PolynomialChaosExpansion : public BasisExpansion {
  friend class PCEFactory;
  
  public:

    PolynomialChaosExpansion(std::shared_ptr<OrthogonalPolynomial>          const& basisCompsIn,
                             std::shared_ptr<muq::Utilities::MultiIndexSet>        multisIn,
                             Eigen::MatrixXd                                const& coeffsIn);

    PolynomialChaosExpansion(std::shared_ptr<OrthogonalPolynomial>          const& basisCompsIn,
                             std::shared_ptr<muq::Utilities::MultiIndexSet>        multisIn,
                             unsigned int                                          outputDim);

    PolynomialChaosExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                             std::shared_ptr<muq::Utilities::MultiIndexSet>            multisIn,
                             Eigen::MatrixXd                                    const& coeffsIn);

    PolynomialChaosExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                             std::shared_ptr<muq::Utilities::MultiIndexSet>            multisIn,
                             unsigned int                                              outputDim);


    virtual ~PolynomialChaosExpansion() = default;


    // ///Print the pce to filename using the << operator.
    // void  Print(std::string filename);
    //
    // ///Print the pce with normalized coefficients to basename+"_pce.dat" using the << operator.
    // void  PrintNormalized(std::string basename);

    ///compute the variance of the current expansion
    Eigen::VectorXd Variance() const;
    Eigen::MatrixXd Covariance() const;
    Eigen::VectorXd Mean() const;

    // /// get the derivative of the variances wrt the polynomial coefficients -- useful in constrained regression
    // Eigen::MatrixXd ComputeVarianceJacobian() const;
    //
    // /** Compute the gradient of the variance in dimension outInd to the polynomial coefficients. */
    // Eigen::VectorXd ComputeVarianceGradient(int outInd) const;
    //
    // /** Compute the Hessian of the variance of a single output dimension with respect to the polynomial coefficients.
    //  *   Since the variance expression is qudratic, this Hessian is constant and thus does not depend on the polynomial
    //  *  coefficients.  For that reason there is no dimension input, even though the output matrix is the Hessian of a
    //  *  single output variance to the polynomial coefficients.
    //  */
    // Eigen::MatrixXd ComputeVarianceHessian() const;


    ///Compute the L2 norm of each output.
    Eigen::VectorXd Magnitude() const;

    /**
     * Compute the weighted sum of polynomial expansions. Slow, because it assumes you haven't been
     * tracking which polynomials to keep, so it does that, then calls the private method. If you
     * know already, your class should be a friend and use the private method directly, as that
     * will be faster.
     * */
    static std::shared_ptr<PolynomialChaosExpansion> ComputeWeightedSum(std::vector<std::shared_ptr<PolynomialChaosExpansion>> expansions,
                                                                        Eigen::VectorXd                                 const& weights);


    ///Compute the Sobol total sensitivity index for the input dimension, for each output dimension
    Eigen::VectorXd TotalSensitivity(unsigned const targetDim) const;

    ///Compute all Sobol total sensitivities. Rows are outputs, each column is an input.
    Eigen::MatrixXd TotalSensitivity() const;

    ///Compute the main sensitivity index for the input dimension, for each output dimension
    Eigen::VectorXd MainSensitivity(unsigned const targetDim) const;

    ///Compute all the main sensitivities. Rows are outputs, each column is an input.
    Eigen::MatrixXd MainSensitivity() const;

    // ///Load an expansion using boost::serialization and return the result
    // static std::shared_ptr<PolynomialChaosExpansion> LoadFromFile(std::string fileName);
    //
    // ///Save an expansion to a file using boost::serialization
    // static void SaveToFile(std::shared_ptr<PolynomialChaosExpansion> expansion, std::string fileName);

  private:

    // PolynomialChaosExpansion(std::vector<std::shared_ptr<OrthogonalPolynomial>> polys,
    //                          unsigned int                                       outputSize);

    /**
     * This function returns the sqrt of the normalization, sqrt(<p*p>),
     * for each PCE basis function p.
     * NB: This function does not change depending on the coefficients.
     */
    Eigen::VectorXd GetNormalizationVec() const;

    /**
     * An internal function to compute the weighted sum of polynomial expansions if you already know which polynomials are
     * included, where polynomials
     * is a clean copy not used by other expansions. Actually does the addition.
     * */
    static std::shared_ptr<PolynomialChaosExpansion> ComputeWeightedSum(std::vector<std::shared_ptr<PolynomialChaosExpansion>> expansions,
                                                                        Eigen::VectorXd const&                                 weights,
                                                                        std::shared_ptr<muq::Utilities::MultiIndexSet> const&  polynomials,
                                                                        std::vector<std::vector<unsigned int>>         const&  locToGlob);

  };

} // namespace muq
} // namespace Approximation



#endif /* POLYNOMIALCHAOSEXPANSION_H_ */
