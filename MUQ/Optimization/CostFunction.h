#ifndef COSTFUNCTION_H_
#define COSTFUNCTION_H_

#include "MUQ/Utilities/VariadicMacros.h"

#include "MUQ/Modeling/ModPiece.h"

namespace muq {
  namespace Optimization {
    class CostFunction : public muq::Modeling::ModPiece {
    public:

      CostFunction(Eigen::VectorXi const& inputSizes);

      virtual ~CostFunction();

      double Cost(muq::Modeling::ref_vector<Eigen::VectorXd> const& input);
      
      template<typename... Args>
	inline double Cost(Args const&... args) {
	return Evaluate(args...).at(0) (0);
      }

      Eigen::VectorXd const& Gradient(unsigned int const inputDimWrt, std::vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity);
      
      template<typename... Args>
	inline Eigen::VectorXd const& Gradient(unsigned int const inputDimWrt, Args const&... args) {
	return ModPiece::Gradient(0, inputDimWrt, args...);
      }
      
    private:

      virtual double CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) = 0;

      virtual void EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

      virtual void GradientImpl(unsigned int const outputDimWrt, unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) override;

      virtual void GradientImpl(unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity);
    };
  } // namespace Optimization
} // namespace muq

#endif
