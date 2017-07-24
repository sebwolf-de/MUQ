#ifndef LINEAROPERATOR_H
#define LINEAROPERATOR_H

#include <Eigen/Core>

#include <memory>
#include <iostream>
#include <exception>

namespace muq
{
namespace Utilities
{


class LinearOperatorTypeException: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Tried creating a linear operator from an unsupported type.  Make sure all necessary headers are included and a child of LinearOperator exists for this type.";
  }
};
    
class LinearOperator;
    
template<typename MatrixType>
struct LinearOperatorFactory
{
    static std::shared_ptr<LinearOperator> Create(MatrixType const& A)
    {
        throw LinearOperatorTypeException();
        
        return std::shared_ptr<LinearOperator>();
        
    };
};

/** @class LinearOperator
 *  @ingroup LinearOperators
 *  @brief Generic linear operator base class
 *  @details In many situations, it is convenient to work with general linear operators instead of specific matrices.
 *      This class provides an abstract base class to help support working with general operators.  

        Children of this class need to implement the Apply and ApplyTranspose pure virtual functions.  
        In addition, each child of this class should provide a template specialization of the static Create
        function (see EigenLinearOperator for an example).

        Typically, an instance of LinearOperator will be created by calling the Create function.  For example:
@code{.cpp}

#include "MUQ/Utilities/LinearAlgebra/EigenLinearOperator.h"

// ...

Eigen::MatrixXd A = Eigen::MatrixXd::Random(10,10);
std::shared_ptr<muq::Utilities::LinearOperator> Aop = muq::Utilities::LinearOperator::Create(A);

@endcode
 */
class LinearOperator {
public:

  LinearOperator(int rowsIn, int colsIn) : ncols(colsIn), nrows(rowsIn) {};

  virtual ~LinearOperator(){};
  
  /** Apply the linear operator to a vector */
  virtual Eigen::MatrixXd Apply(Eigen::Ref<const Eigen::MatrixXd> const& x) = 0;
  
  /** Apply the transpose of the linear operator to a vector. */
  virtual Eigen::MatrixXd ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x) = 0;

  /** Fills in the reference \f$y\f$ with \f$y=Ax\f$ */
  virtual void Apply(Eigen::Ref<const Eigen::MatrixXd> const& x, Eigen::Ref<Eigen::MatrixXd> y);

  /** Fill in the reference \f$y\f$ with \f$y = A^Txf$ */
  virtual void ApplyTranspose(Eigen::Ref<const Eigen::MatrixXd> const& x, Eigen::Ref<Eigen::MatrixXd> y);
  
  /** The output dimension of the linear operator. */
  int rows() const { return nrows; }

  /** The input dimension of the linear operator. */
  int cols() const { return ncols; }

  /** Get a dense matrix representation of this linear operator. */
  virtual Eigen::MatrixXd GetMatrix();
  

  template<typename OtherType>
  static std::shared_ptr<LinearOperator> Create(OtherType const& A)
  {
      return LinearOperatorFactory<OtherType>::Create(A);
  }
  
protected:

  const int ncols;
  const int nrows;
};




} // namespace Utilities
} // namespace MUQ

#endif
