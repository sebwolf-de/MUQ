#ifndef COVARIANCEKERNELS_H_
#define COVARIANCEKERNELS_H_

#include "MUQ/Approximation/GaussianProcesses/ConcatenateKernel.h"
#include "MUQ/Approximation/GaussianProcesses/ConstantKernel.h"
#include "MUQ/Approximation/GaussianProcesses/CoregionalKernel.h"
#include "MUQ/Approximation/GaussianProcesses/LinearTransformKernel.h"
#include "MUQ/Approximation/GaussianProcesses/PeriodicKernel.h"
#include "MUQ/Approximation/GaussianProcesses/ProductKernel.h"
#include "MUQ/Approximation/GaussianProcesses/SquaredExpKernel.h"
#include "MUQ/Approximation/GaussianProcesses/SumKernel.h"
#include "MUQ/Approximation/GaussianProcesses/WhiteNoiseKernel.h"


namespace muq
{
namespace Approximation
{
    
/**
\defgroup CovarianceKernels Covariance Kernels
\ingroup GaussianProcesses

Gaussian processes are defined by their mean function and covariance kernel.  In this group, we provide tools for constructing and using covariance kernels.  Several simple kernels are provided as well as tools for combining simple kernels into more complicated kernels.

C++11 introduced the "auto" keyword that tells the compiler to try to automatically detect the type of a variable.  For example, consider a kernel of the form
\f[
k(x,y) = k_1(x,y) * k_2(x,y) + k_3(x,y),
\f]
where \f$k_1\f$ is a squared exponential kernel, \f$k_2\f$ is a periodic kernel and \f$k_3\f$ is a white noise kernel.  In code, this can be implemented as
\code{.cpp}
auto kernel1 = SquaredExpKernel(var1, length1);
auto kernel2 = PeriodicKernel(var2, length2, period);
auto kernel3 = WhiteNoiseKernel(var3);

auto kernel = kernel1*kernel2 + kernel3;
\endcode
or, more succinctly, as 
\code{.cpp}
auto kernel = SquaredExpKernel(var1, length1) * PeriodicKernel(var2, length2, period) + WhiteNoiseKernel(var3);
\endcode

In either case, the kernel variable is an instance of "SumKernel<ProductKernel<SquaredExpKernel,PeriodicKernel>, WhiteNoiseKernel>".  The "auto" keyword allows us to avoid typing this long type.
 
*/


} // namespace Approximation
} // namespace muq



#endif // #ifndef COVARIANCEKERNELS_H_
