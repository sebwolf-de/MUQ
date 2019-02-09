#include "MUQ/Approximation/Quadrature/ClenshawCurtisQuadrature.h"

using namespace muq::Approximation;


ClenshawCurtisQuadrature::ClenshawCurtisQuadrature() : Quadrature(1) {}


void ClenshawCurtisQuadrature::Compute(unsigned int index) {

  unsigned int order  = index+1;

  pts.resize(1,order);
  wts.resize(order);

  if(order==1){
    pts(0,0) = 0.0;
    wts(0) = 2.0;
    return;
  }

  for(int i=0; i<order; ++i)
    pts(0,i) = std::cos( double(order-i-1) * pi / double(order - 1));

  pts(0,0) = -1.0;
  if((order%2)==1)
    pts((order+1)/2-1) = 0.0;
  pts(order-1) = 1.0;

  wts = Eigen::VectorXd::Ones(order);

  for(int i=0; i<order; ++i) {
    double theta = double(i) * pi / double(order-1);

    for(int j=0; j<(order-1)/2; ++j){

      double b;
      if(2*(j+1)==(order-1)){
        b = 1.0;
      }else{
        b = 2.0;
      }

      wts(i) -= b*std::cos(2.0*(j+1)*theta) / (4.0*std::pow(j+1,2) - 1);
    }
  }

  wts(0) /= double(order-1);
  wts.segment(1,order-2) *= 2.0/(order-1);
  wts(order-1) /= double(order-1);

}
