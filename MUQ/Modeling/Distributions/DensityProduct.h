#ifndef DENSITYPRODUCT_H
#define DENSITYPRODUCT_H

#include "MUQ/Modeling/Distributions/Density.h"

namespace muq {
  namespace Modeling {

    class DensityProduct : public DensityBase {
    public:

      /** Constructs a density product accepting numPiecesIn inputs.  For example,
          a Bayesian posterior density is composed of a likelihood and prior, so
          numPiecesIn=2.
          @param[in] numPiecesIn Number of input densities.
      */
      DensityProduct(int numPiecesIn);

      virtual ~DensityProduct() = default;

    private:

      virtual double LogDensityImpl(ref_vector<boost::any> const& inputs) override;

      virtual void JacobianImpl(unsigned int           const  wrtIn,
                                unsigned int           const  wrtOut,
                                ref_vector<boost::any> const& inputs) override;

      virtual void JacobianActionImpl(unsigned int           const  wrtIn,
                                      unsigned int           const  wrtOut,
                                      boost::any             const& vec,
                                      ref_vector<boost::any> const& inputs) override;

      virtual void JacobianTransposeActionImpl(unsigned int           const  wrtIn,
                                               unsigned int           const  wrtOut,
                                               boost::any             const& vec,
                                               ref_vector<boost::any> const& inputs) override;
                                               
    };
  } // namespace Modeling
} // namespace muq

#endif
