#include "MUQ/Modeling/Distributions/Uniform.h"

#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

Uniform::~Uniform() {}

Uniform::Uniform(std::vector<std::pair<double, double> > const& bounds) : Distribution(), bounds(ComputeBounds(bounds)) {
  // initialize the any algebra
  algebra = std::make_shared<AnyAlgebra>();
}

std::vector<std::pair<double, double> > Uniform::CreateBounds(std::vector<std::pair<double, double> >& bounds, std::pair<double, double> last) {
  // add the last pair
  bounds.push_back(last);

  // return the upper/lower bound paris
  return bounds;
}

std::vector<std::pair<double, double> > Uniform::ComputeBounds(std::vector<std::pair<double, double> > const& bounds) {
  std::vector<std::pair<double, double> > new_bounds(bounds);
  
  for( auto it=new_bounds.begin(); it!=new_bounds.end(); ++it ) {
    it->second = it->second-it->first;
  }

  return new_bounds;
}

double Uniform::LogDensityImpl(ref_vector<boost::any> const& inputs) const {
  // get the dimension
  const unsigned int dim = algebra->Size(inputs[0].get());
  assert(dim==bounds.size());

  // loop though each dimension
  for( unsigned int i=0; i<dim; ++i ) {
    // get the value of the point for this dimension
    const double x = boost::any_cast<double const>(algebra->AccessElement(inputs[0].get(), i));

    // if this dimension is outside the hypercube
    if( x<bounds[i].first || x>bounds[i].first+bounds[i].second ) {
      // return negative infinity
      return -std::numeric_limits<double>::infinity();
    }
  }

  // the point is inside the domain...
  return 1.0;
}

boost::any Uniform::SampleImpl(ref_vector<boost::any> const& inputs) {
  assert(bounds.size()>0);
  
  if( bounds.size()==1 ) { // one dimensional
    // sample and return a double
    return bounds[0].first+bounds[0].second*RandomGenerator::GetUniform();
  }

  // create a vector
  Eigen::VectorXd sample(bounds.size());
  for( unsigned int i=0; i<bounds.size(); ++i ) {
    // sample each dimension
    sample(i) = bounds[i].first+bounds[i].second*RandomGenerator::GetUniform();
  }
  
  return sample;
}
