#ifndef LOBPCG_H
#define LOBPCG_H


#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

#include <boost/property_tree/ptree_fwd.hpp>
#include <Eigen/Core>

namespace muq{
namespace Modeling{

  /** @class LOBPCG
      @ingroup LinearOperators
      @brief The Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG) method for matrix-free computation of eigenvalues and eigenvectors.
      @details This class solves generalized eigenvalue problems of the form \f$Av = \lambda Bv\f$ when the matrix \f$A\f$ is symmetric and the matrix \f$Bf$ is symmetric positive definite.  It uses the The Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG) method described in
      "TOWARD THE OPTIMAL PRECONDITIONED EIGENSOLVER: LOCALLY OPTIMAL BLOCK PRECONDITIONED CONJUGATE GRADIENT METHOD" by ANDREW V. KNYAZEV.
   */
  class LOBPCG
  {
  public:

    /**
    @param[in] numEigsIn The maximum number of eigenvalues to compute.
    @param[in] tolIn Solver tolerance
    @param[in] maxItsIn The maximum number of iterations the solver is allowed to take.
    @param[in] largestIn A boolean specifying if the largest (true) or smallest (false) eigenvalues should be computed.
    @param[in] verbosityIn An integer {0,1,2,3} specifying how much information the solver should print.  When verbosityIn=0 (default), nothing is printed.
    */
    LOBPCG(int    numEigsIn,
           double tolIn=-1,
           int    maxItsIn=-1,
           bool   largestIn=true,
           int    verbosityIn=0);

    /**
      <B>Configuration Parameters:</B>
      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "NumEigs"     | integer       | -             | The maximum number of eigenvalues to compute. |
      "Tolerance"   | double        | \f$d\sqrt{\epsilon}\f$ | Tolerance used for stopping criteria.  |
      "MaxIts"      | integer       | min(d,20)   | Maximum number of iterations to take. |
      "Largest"     | bool          | True |  If true, the largest eigenvalues are computed.  If false, the smallest eigenvalues are returned. |
      "Verbosity"   | integer       | 0  |  Specifies how much is printed to the screen.  Valid values are 0 (print nothing), 1, 2, or 3       |
    */
    LOBPCG(boost::property_tree::ptree const& options);

    /**
    Compute the generalized eigenvalues and eigenvectors of a symmetric system \f$Av = \lambda Bv\f$.  If compute has been previously called, the eigenvalues and eigenvectors from the previous call will be used as an initial guess for this call to the solver.  This can speed up the solver if A changes slightly and the eigenvalues need to be recomputed (as in DILI MCMC).

    @param[in] A A shared pointer to a muq::Modeling::LinearOperator instance specifying the matrix A.
    @param[in] B A shared pointer to a muq::Modeling::LinearOperator instance specifying the matrix B.  If not specified or set to nullptr, B will be set to the identity.
    @param[in] M An optional preconditioner.  If specified, the LinearOperator M should approximate the inverse of A.
    */
    LOBPCG& compute(std::shared_ptr<LinearOperator> const& A,
                    std::shared_ptr<LinearOperator>        B = nullptr,
                    std::shared_ptr<LinearOperator>        M = nullptr);

    /**
    Compute the generalized eigenvalues and eigenvectors of a symmetric system \f$Av = \lambda Bv\f$ that are in the B-orthogonal complement to some constraints \f$Y\f$.  More precisely, the computed eigenvectors \f$v_i\v$ will satisfy \f$Y^T B v_i=0\f$.
    @param[in] A A shared pointer to a muq::Modeling::LinearOperator instance specifying the matrix A.
    @param[in] constMat A matrix defining the constraints.  This matrix should have the same number of rows as A.
    @param[in] B A shared pointer to a muq::Modeling::LinearOperator instance specifying the matrix B.  If not specified or set to nullptr, B will be set to the identity.
    @param[in] M An optional preconditioner.  If specified, the LinearOperator M should approximate the inverse of A.
    */
    LOBPCG& compute(std::shared_ptr<LinearOperator> const& A,
                    Eigen::MatrixXd                 const& constMat,
                    std::shared_ptr<LinearOperator>        B = nullptr,
                    std::shared_ptr<LinearOperator>        M = nullptr);

    /** Return a reference to the computed vector of eigenvalues.  The vector
        will only be valid after calling compute.
    */
    Eigen::VectorXd const& eigenvalues() const{return eigVals;}

    /** Return a matrix whose columns contain the computed eigenvectors.  The
        matrix will only be valid after calling compute.
    */
    Eigen::MatrixXd const& eigenvectors() const{return eigVecs;};

    /** Resets the current eigenvalues and eigenvectors so that a random initial
        guess is used during the next call to compute instead of reusing the
        previously computed eigenvalues and eigenvectors as initial guesses.
        @param[in] dim The dimension of the problem
     */
    LOBPCG& reset(int dim);

  private:

    /** Makes the columns of a matrix V orthonormal wrt the B inner product \f$v^T B v\f$. */
    class Orthonormalizer{
    public:
      Orthonormalizer(std::shared_ptr<LinearOperator> const& Bin) : B(Bin){};

      void ComputeInPlace(Eigen::MatrixXd& V);
      void ComputeInPlace(Eigen::MatrixXd& V, Eigen::MatrixXd const& BVin);

      Eigen::MatrixXd Compute(Eigen::MatrixXd const& V);
      Eigen::MatrixXd Compute(Eigen::MatrixXd const& V, Eigen::MatrixXd const& BVin);

      Eigen::MatrixXd InverseVBV() const;
      Eigen::MatrixXd const& GetBV() const{return BV;};
      Eigen::MatrixXd& GetBV(){return BV;};

      int vDim;
      std::shared_ptr<LinearOperator> B;
      Eigen::MatrixXd BV;
      Eigen::MatrixXd VBV_Chol;
    };

    class Constraints{
    public:
      Constraints(std::shared_ptr<LinearOperator> const& B,
                  Eigen::MatrixXd                 const& constVec);

      void ApplyInPlace(Eigen::MatrixXd& x);

      int size() const{return BY.cols();};

    private:
      Eigen::MatrixXd BY;
      Eigen::MatrixXd const& Y;
      Eigen::LLT<Eigen::MatrixXd> YBY_llt;

    }; // class LOBPCG::Constraints

    /// The maximum number of eigenvalues to compute
    int numEigs;

    /// Solver tolerance
    double tol;

    /// Maximum number of iterations
    int maxIts;

    /// Do we want to pursue the largest or smallest eigenvalues?  If largest==true, then we go for the big ones.
    const bool largest;

    /// Controls how much information we want to print
    int verbosity;

    Eigen::VectorXd eigVals;
    Eigen::MatrixXd eigVecs;

    std::shared_ptr<LinearOperator> A, B, M;

  }; // class LOBPCG
}
}



#endif // LOBPCG_H
