#include <MUQ/Modeling/Distributions/Distribution.h>

class Target : public muq::Modeling::Distribution {
public:

  inline Target(double const p) : muq::Modeling::Distribution(1), p(p) {}

  inline virtual ~Target() {}

  inline virtual double LogDensityImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs) override {
    // get the point where we are evaluating the log density
    const double x = inputs.at(0) (0);

    return -std::pow(std::fabs(x), p)+std::sin(2.0*M_PI*x);
  }

private:

  /// The power that determines tail behavior
  const double p;
};
