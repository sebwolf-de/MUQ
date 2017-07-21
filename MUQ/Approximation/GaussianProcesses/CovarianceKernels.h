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
#include "MUQ/Approximation/GaussianProcesses/MaternKernel.h"

namespace muq
{
namespace Approximation
{
    
/**
\defgroup CovarianceKernels Covariance Kernels
\ingroup GaussianProcesses
\brief Covariance kernels for defining Gaussian Processes

Gaussian processes are defined by their mean function and covariance kernel.  In this group, we provide tools for constructing and using covariance kernels.  Several simple kernels are provided as well as tools for combining simple kernels into more complicated kernels.  

Templates are used extensively in this module to maximize performance (by reducing virtual function calls).  Thus, using the c++11 "auto" keyword can result in much cleaner and readable code.

<h2>Combining Kernels</h2>
As an example of constructing a covariance kernel that is constructed from several simple kernels, consider a kernel of the form
\f[
k(x_1,x_2) = k_1(x_1,x_2)\, k_2(x_1,x_2) + k_3(x_1,x_2),
\f]
where \f$k_1\f$ is a squared exponential kernel, \f$k_2\f$ is a periodic kernel and \f$k_3\f$ is a white noise kernel.  Assume the dimension of \f$x_1\f$ and \f$x_2\f$ is stored in a variable called <code>dim</code>.  Then, the combined kernel \f$k(x,y)\f$ can be constructed with the following snippet:
\code{.cpp}
#include "MUQ/Approximation/GaussianProcesses/CovarianceKernels.h"

// Other setup ...

auto k1 = SquaredExpKernel(dim, var1, length1);
auto k2 =   PeriodicKernel(dim, var2, length2, period);
auto k3 = WhiteNoiseKernel(dim, var3);

auto k = kernel1*kernel2 + kernel3;
\endcode
or, more succinctly, as 
\code{.cpp}
auto k = SquaredExpKernel(dim, var1, length1) * PeriodicKernel(dim, var2, length2, period) + WhiteNoiseKernel(dim, var3);
\endcode

In either case, the <code>k</code> variable is an instance of "SumKernel<ProductKernel<SquaredExpKernel,PeriodicKernel>, WhiteNoiseKernel>".  The "auto" keyword allows us to avoid typing this long type.
 
<h2>Anisotropic Kernels</h2>
In many cases, the correlation of a GP in one dimension will be different that the correlation is some other direction.  For example, let's say we have two spatial variables \f$x\f$ and \f$y\f$ as well as a kernel of the form
\f[
k([x_1,y_1], [x_2, y_2]) = k_x(x_1, x_2)\, k_y(y_1, y_2). 
\f]
Such kernels commonly arise when modeling anisotropic media (e.g., hydraulic conductivity fields).  In MUQ, it is possible to specify the dimensions that are used by a kernel.  For example, if \f$k_1\f$ and \f$k_2\f$ were both squared exponential kernels, than \f$k\f$ could be defined as
\code{.cpp}

// Only keep the 0 index
std::vector<unsigned> indsx = {0};
auto kx = SquaredExpKernel(2, indsx, varx, Lx);

// Only keep the 1 index
std::vector<unsigned> inds2 = {1};
auto ky = SquaredExpKernel(2, indsy, vary, Ly);

auto k = k1 * k2;

\endcode

It is also possible to define more complicated relationships.  For example, consider a third component \f$z\f$, and let the kernel be defined as 
\f[
k([x_1,y_1,z_1], [x_2, y_2, z_2]) = k_{xy}([x_1,y_1], [x_2,y_2])\, k_z(z_1, z_2). 
\f]
This kernel might be constructed with
\code{.cpp}

// Keep both the 0 and 1 indices
std::vector<unsigned> indsxy = {0,1};
auto kx = SquaredExpKernel(3, indsxy, varxy, Lxy);

// Only keep the 2 index
std::vector<unsigned> inds2 = {2};
auto kz = SquaredExpKernel(3, indsz, varz, Lz);

auto k = kxy * kz;

\endcode

<h2>Vector-Valued Kernels</h2>
Gaussian processes are used most often to characterize scalar-valued functions.  However, gaussian processes can also be incredibly useful for characterizing vector-valued functions.  A common way to handle vector-valued Gaussian processes, and the one employed by MUQ, is to adopt a ``Linear Model of Coregionalization."  In this setting, each component of the vector-valued Gaussian process is expressed as a sum of independent latent Gaussian processes. See \cite Wackernagel2010 for more details.

In MUQ, vector-valued Gaussian processes are represented through a vector-valued mean function and a matrix-valued covariance function.

A common way of constructing a vector-valued Gaussian process is to first define a \f$D\times D\f$ dimensional marginal covariance \f$\Sigma\f$ and define the process with scalar-valued Gaussian processes on the eigenvectors of \f$\Sigma\f$.  Let \f$Q\Lambda Q^T=\Sigma\f$ be the eigenvalue decomposition of the covariance matrix, where \f$Q\f$ is a matrix whose columns contain the eigenvectors and \f$\Lambda\f$ is a diagonal matrix containing the corresponding eigenalues.  With this decomposition, the muq::Approximation::CoregionalKernel defines the matrix-valued kernel \f$k\f$ as 
\f[
k(x_1,x_2) = \sum_{i=1}^D \sqrt{\lambda_i} k_i(x_1,x_2) Q_i Q_i^T,
\f]
where each \f$k_i\f$ is a scalar-valued kernel.  

In MUQ, such a kernel can be constructed using
\code{.cpp}

Eigen::MatrixXd marginalCov(2,2);
marginalCov << 1.0, 0.8,
               0.8, 1.0;

// The kernel along the first eigenvector of marginalCov
auto k1 = SquaredExpKernel(2, var1, length1);

// The kernel along the second eigenvector of marginalCov
auto k2 = SquaredExpKernel(2, var1, length1);

auto k = CoregionTie(marginalCov, k1, k2);

\endcode
Note that we use Eigen to compute the eigendecomposition of the covariance.  Eigen sorts eigenvalues in ascending order, meaning that the kernel <code>k1</code> in this snippet will correpsond to te smallest eigenvalue.

The constructor of the muq::Approximation::CoregionalKernel can also be called directly.  For example, the same kernel can be constructed using
\code{.cpp}

Eigen::MatrixXd marginalCov(2,2);
marginalCov << 1.0, 0.8,
               0.8, 1.0;

std::vector<std::shared_ptr<KernelBase>> ks(2);
ks.at(0) = std::make_shared<SquaredExpKernel>(2, var1, length1);
ks.at(1) = std::make_shared<SquaredExpKernel>(2, var2, length2);

CoregionalKernel k(2, marginalCov, kernels);

\endcode


<h3> Linear Operations </h3>
Suppose we have a vector-valued Gaussian Process \f$f\sim GP( \mu(x), k_f(x, x^\prime) )\f$ where \f$\mu(x)\f$ returns a vector with \f$N\f$ components and \f$k(x,x^\prime)\f$ returns an \f$N\times N\f$ marginal covariance matrix.  Now, let \f$A\f$ be an \f$M\times N\f$ matrix and consider the new Gaussian process \f$g = Af\f$.  The mean function of \f$g\f$ is given by \f$A \mu(x)\f$ and the covariance kernel is given by \f$k_g(x,x^\prime) = A k_f(x, x^\prime) A^T \f$.  Constructing this type of kernel with MUQ is straightforward:
\code{.cpp}

// First set up a vector-valued kernel (a coregional Kernel in this case)
Eigen::MatrixXd marginalCov(2,2);
marginalCov << 1.0, 0.8,
               0.8, 1.0;

auto k1 = SquaredExpKernel(2, var1, length1);
auto k2 = SquaredExpKernel(2, var1, length1);

auto kf = CoregionTie(marginalCov, k1, k2);

// Now apply the linear operation
Eigen::MatrixXd A(3,2);
A << 1.0, 1.0,
     0.5, 0.75,
     1.0, 0.0;

auto kg = A * kf;

\endcode

*/


} // namespace Approximation
} // namespace muq



#endif // #ifndef COVARIANCEKERNELS_H_
