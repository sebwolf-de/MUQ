#include "MUQ/Modeling/Distributions/DensityProduct.h"
#include "MUQ/Utilities/AnyHelpers.h"


using namespace muq::Modeling;
using namespace muq::Utilities;

DensityProduct::DensityProduct(int numPiecesIn) : DensityBase()
{
  numInputs = numPiecesIn;
}


double DensityProduct::LogDensityImpl(ref_vector<boost::any> const& inputs)
{
  double sum = 0.0;
  for(int i=0; i<inputs.size(); ++i){
    double part = AnyConstCast(inputs.at(i).get());
    sum += part;
  }
  return sum;
}

void DensityProduct::JacobianImpl(unsigned int           const  wrtIn,
                                  unsigned int           const  wrtOut,
                                  ref_vector<boost::any> const& inputs)
{
  jacobian = Eigen::MatrixXd::Ones(1,1).eval();
};

void DensityProduct::JacobianActionImpl(unsigned int           const  wrtIn,
                                        unsigned int           const  wrtOut,
                                        boost::any             const& vec,
                                        ref_vector<boost::any> const& inputs)
{
  jacobianAction = vec;
};

void DensityProduct::JacobianTransposeActionImpl(unsigned int           const  wrtIn,
                                                 unsigned int           const  wrtOut,
                                                 boost::any             const& vec,
                                                 ref_vector<boost::any> const& inputs)
{
  jacobianTransposeAction = vec;
};
