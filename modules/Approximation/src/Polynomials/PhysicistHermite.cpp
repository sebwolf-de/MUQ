#include "MUQ/Approximation/Polynomials/PhysicistHermite.h"

using namespace muq::Approximation;

PhysicistHermite::PhysicistHermite() : Polynomial() {}

PhysicistHermite::~PhysicistHermite() {}

double PhysicistHermite::alpha(unsigned int k, double x) const {
  return -2.0*x;
}

double PhysicistHermite::beta(unsigned int k, double x) const {
  return 2.0*k;
}

double PhysicistHermite::phi0(double x) const {
  return 1.0;
}

double PhysicistHermite::phi1(double x) const {
  return 2.0*x;
}

REGISTER_POLYNOMIAL_FAMILY(PhysicistHermite)
