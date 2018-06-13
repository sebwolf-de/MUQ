#ifndef PYMODPIECE_H
#define PYMODPIECE_H

#include "MUQ/Modeling/ModPiece.h"

namespace muq {
  namespace Modeling {

    class PyModPiece : public ModPiece {

    public:
      PyModPiece(Eigen::VectorXi const& inputSizes,
                 Eigen::VectorXi const& outputSizes);

    protected:
      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) override;

      virtual void EvaluateImpl(std::vector<Eigen::VectorXd> const& input) = 0;

      virtual void GradientImpl(unsigned int                const  outputDimWrt,
                                unsigned int                const  inputDimWrt,
                                ref_vector<Eigen::VectorXd> const& input,
                                Eigen::VectorXd             const& sensitivity) override;

      virtual void GradientImpl(unsigned int                 const  outputDimWrt,
                                unsigned int                 const  inputDimWrt,
                                std::vector<Eigen::VectorXd> const& input,
                                Eigen::VectorXd              const& sensitivity) = 0;

      virtual void JacobianImpl(unsigned int                const  outputDimWrt,
                                unsigned int                const  inputDimWrt,
                                ref_vector<Eigen::VectorXd> const& input) override;

      virtual void JacobianImpl(unsigned int                 const  outputDimWrt,
                                unsigned int                 const  inputDimWrt,
                                std::vector<Eigen::VectorXd> const& input) = 0;

      virtual void ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                     unsigned int                const  inputDimWrt,
                                     ref_vector<Eigen::VectorXd> const& input,
                                     Eigen::VectorXd             const& vec) override;

      virtual void ApplyJacobianImpl(unsigned int                 const  outputDimWrt,
                                     unsigned int                 const  inputDimWrt,
                                     std::vector<Eigen::VectorXd> const& input,
                                     Eigen::VectorXd              const& vec) = 0;

    private:
      std::vector<Eigen::VectorXd> ToStdVec(ref_vector<Eigen::VectorXd> const& input);

    };
  }
}

#endif
