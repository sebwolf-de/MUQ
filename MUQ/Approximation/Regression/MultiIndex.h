#ifndef MULTIINDEX_H_
#define MULTIINDEX_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Approximation {
    /// A multi-index determines the powers of a multi-dimensional polynomial
    /**
       Let the \f$n\f$-dimensional multi-index be defined by
       \f{eqnarray*}{
       \alpha = \left[\alpha_1,\, ...,\, \alpha_n\right]
       \f}
       Such that the magnitude is 
       \f{eqnarray*}{
       \vert\alpha\vert = \sum_{i=1}^n \alpha_i
       \f}
       and the power is
       \f{eqnarray*}{
       x^{\alpha} = \prod_{i=1}^n p_{\alpha_i}(x_i)
       \f}
       where \f$p_{\alpha_i}(x_i)\f$ is a \f$\alpha_i^{th}\f$ degree polynomial (e.g., a monomial, Hermite, or Legendre polynomial).
     */
    class MultiIndex : public muq::Modeling::WorkPiece {
    public:

      /**
	 Creates a set of all multi-indexes \f$\alpha\f$ whose order is \f$L_l \leq \vert\alpha\vert \leq L_u\f$ (\f$L_l\f$ and \f$L_u\f$ are the lower and upper bounds of the multi-index order).

	 @param[in] n The dimension of the multi-index
	 @param[in] maxOrder The maximum order of the multi-index
	 @param[in] minOrder The minimum order of the multi-index (defaults to zero)
       */
      MultiIndex(unsigned int const n, unsigned int const maxOrder, unsigned int const minOrder = 0);

      virtual ~MultiIndex();

      /// Return the number of mult-indices in this set
      /**
	 \return The size of muq::Approximation::MultiIndex::alphas
       */
      unsigned int Size() const;
      
    private:

      /**
	 @param[in] alpha The multi index that we may or may not add to the set of multi-indices
	 @param[in] minOrder The minimum order of the multi-index (defaults to zero)
	 @param[in] maxOrder The maximum order of the multi-index
	 @param[in] dim We have already added all of the multi-indices up to this dimension (defaults to zero)
      */
      void AddMulti(Eigen::VectorXi const alpha, unsigned int const minOrder, unsigned int const maxOrder, unsigned int const dim = 0);

      /**
	 inputs:
	 <ol>
	 <li> The \f$i^{th}\f$ index---we want to retrive this multi-index (unsigned int)
	 </ol>
	 outputs:
	 <ol>
	 <li> The \f$i^{th}\f$ multi-index (Eigen::VectorXi)
	 </ol>
       */
      void EvaluateImpl(muq::Modeling::ref_vector<boost::any> const& inputs) override;

      /// A list of all of the mutli-indices with the appropriate order
      std::vector<Eigen::VectorXi> alphas;
      
    };
  } // namespace Approximation
} // namespace muq

#endif
