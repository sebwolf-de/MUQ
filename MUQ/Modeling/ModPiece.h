#ifndef MODPIECE_H
#define MODPIECE_H

#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Utilities/VariadicMacros.h"

#include <Eigen/Core>
#include <vector>

namespace muq{
  namespace Modeling{

  class ModPiece : public WorkPiece
  {

  public:

    ModPiece(Eigen::VectorXi const& inputSizes,
             Eigen::VectorXi const& outputSizes);

    virtual ~ModPiece() = default;


    /** @brief Get the average run time for one of the implemented methods.
     *   @details This function returns the average wall clock time (in milliseconds) for the EvaluateImpl, GradientImpl,
     * JacobianImpl, JacobianActionImpl, or HessianImp functions.
     *            If the function was never called, -1 is returned.
     *   @param[in] method The implemented function of interest.  Possible options are "Evaluate", "Gradient", "Jacobian",
     * "JacobianAction", or "Hessian"
     *   @return The average wall clock time in milli-seconds.
     */
    virtual double GetRunTime(const std::string& method="Evaluate") const override;
    virtual void ResetCallTime() override;

    /** @brief get the number of times one of the implemented methods has been called.
     *   @details This function returns the number of times the EvaluateImpl, GradientImpl, JacobianImpl,
     * JacobianActionImpl, or HessianImp functions have been called.
     *   @param[in] method The implemented function of interest.  Possible options are "Evaluate", "Gradient", "Jacobian",
     * "JacobianAction", or "Hessian"
     *   @return An integer with the number of calls.
     */
    virtual unsigned long int GetNumCalls(const std::string& method = "Evaluate") const override;


    virtual std::vector<Eigen::VectorXd> const& Evaluate(std::vector<Eigen::VectorXd> const& input);
    virtual std::vector<Eigen::VectorXd> const& Evaluate(ref_vector<Eigen::VectorXd> const& input);
    VARIADIC_TO_REFVECTOR(Evaluate, Eigen::VectorXd, std::vector<Eigen::VectorXd> const&);


    virtual Eigen::VectorXd const& Gradient(unsigned int                 const  outputDimWrt,
                                            unsigned int                 const  inputDimWrt,
                                            std::vector<Eigen::VectorXd> const& input,
                                            Eigen::VectorXd              const& sensitivity);

    virtual Eigen::VectorXd const& Gradient(unsigned int                const  outputDimWrt,
                                            unsigned int                const  inputDimWrt,
                                            ref_vector<Eigen::VectorXd> const& input,
                                            Eigen::VectorXd             const& sensitivity);

    inline Eigen::VectorXd const& Gradient(unsigned int outWrt, unsigned int inWrt, Eigen::VectorXd const& last, Eigen::VectorXd const& sens) {     \
      ref_vector<Eigen::VectorXd> vec;
      vec.push_back(std::cref(last));                                                                                               \
      return Gradient(outWrt, inWrt, vec, sens);                                                                              \
    }
    template<typename... Args>
    inline Eigen::VectorXd  const& Gradient(unsigned int wrtOut, unsigned int wrtIn, Args const&... args) {
      ref_vector<Eigen::VectorXd> vec;
      return GradientRecurse(wrtOut, wrtIn, vec, args...);
    }


    virtual Eigen::MatrixXd const& Jacobian(unsigned int                 const  outputDimWrt,
                                            unsigned int                 const  inputDimWrt,
                                            std::vector<Eigen::VectorXd> const& input);

    virtual Eigen::MatrixXd const& Jacobian(unsigned int                const  outputDimWrt,
                                            unsigned int                const  inputDimWrt,
                                            ref_vector<Eigen::VectorXd> const& input);


    template<typename... Args>
    inline Eigen::MatrixXd const& Jacobian(unsigned int outWrt, unsigned int inWrt, Args const&... args) {
      ref_vector<Eigen::VectorXd> vec;
      return Jacobian(outWrt, inWrt, vec, args...);
    }
    template<typename... Args>
    inline Eigen::MatrixXd JacobianByFD(unsigned int outWrt, unsigned int inWrt, Args const&... args) {
      ref_vector<Eigen::VectorXd> vec;
      return JacobianByFD(outWrt, inWrt, vec, args...);
    }
    template<typename... Args>
    inline Eigen::MatrixXd ApplyJacobianByFD(unsigned int outWrt, unsigned int inWrt, Args const&... args) {
      ref_vector<Eigen::VectorXd> vec;
      return ApplyJacobianByFD(outWrt, inWrt, vec, args...);
    }


    virtual Eigen::VectorXd const& ApplyJacobian(unsigned int                 const  outputDimWrt,
                                                 unsigned int                 const  inputDimWrt,
                                                 std::vector<Eigen::VectorXd> const& input,
                                                 Eigen::VectorXd              const& vec);

    virtual Eigen::VectorXd const& ApplyJacobian(unsigned int                const  outputDimWrt,
                                                 unsigned int                const  inputDimWrt,
                                                 ref_vector<Eigen::VectorXd> const& input,
                                                 Eigen::VectorXd             const& vec);

    // virtual Eigen::VectorXd ApplyHessian(int                         const  outputDimWrt,
    //                                      int                         const  inputDimWrt1,
    //                                      int                         const  inputDimWrt2,
    //                                      ref_vector<Eigen::VectorXd> const& input,
    //                                      Eigen::VectorXd             const& sensitivity
    //                                      Eigen::VectorXd             const& vec);



    virtual Eigen::VectorXd GradientByFD(unsigned int                 const  outputDimWrt,
                                         unsigned int                 const  inputDimWrt,
                                         std::vector<Eigen::VectorXd> const& input,
                                         Eigen::VectorXd              const& sensitivity);

    virtual Eigen::VectorXd GradientByFD(unsigned int                const  outputDimWrt,
                                         unsigned int                const  inputDimWrt,
                                         ref_vector<Eigen::VectorXd> const& input,
                                         Eigen::VectorXd             const& sensitivity);

    virtual Eigen::MatrixXd JacobianByFD(unsigned int                 const  outputDimWrt,
                                         unsigned int                 const  inputDimWrt,
                                         std::vector<Eigen::VectorXd> const& input);

    virtual Eigen::MatrixXd JacobianByFD(unsigned int                const  outputDimWrt,
                                         unsigned int                const  inputDimWrt,
                                         ref_vector<Eigen::VectorXd> const& input);

    virtual Eigen::VectorXd ApplyJacobianByFD(unsigned int                 const  outputDimWrt,
                                              unsigned int                 const  inputDimWrt,
                                              std::vector<Eigen::VectorXd> const& input,
                                              Eigen::VectorXd              const& vec);

    virtual Eigen::VectorXd ApplyJacobianByFD(unsigned int                const  outputDimWrt,
                                              unsigned int                const  inputDimWrt,
                                              ref_vector<Eigen::VectorXd> const& input,
                                              Eigen::VectorXd             const& vec);

    const Eigen::VectorXi inputSizes;
    const Eigen::VectorXi outputSizes;

  protected:

    // The following variables keep track of how many times the Implemented functions, i.e. EvaluateImpl, GradientImpl,
    // etc... are called.
    unsigned long int numGradCalls   = 0;
    unsigned long int numJacCalls    = 0;
    unsigned long int numJacActCalls = 0;
    unsigned long int numHessCalls   = 0;

    // these variables keep track of the total wall-clock time spent in each of the Implemented functions.  They are in
    // units of milliseconds
    double gradTime   = 0;
    double jacTime    = 0;
    double jacActTime = 0;
    double hessTime   = 0;

    std::vector<Eigen::VectorXd> outputs;
    Eigen::VectorXd gradient;
    Eigen::VectorXd jacobianAction;
    Eigen::MatrixXd jacobian;

    void CheckInputs(ref_vector<Eigen::VectorXd> const& input);

    virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override;

    virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) = 0;

    virtual void GradientImpl(unsigned int                const  outputDimWrt,
                              unsigned int                const  inputDimWrt,
                              ref_vector<Eigen::VectorXd> const& input,
                              Eigen::VectorXd             const& sensitivity);

    virtual void JacobianImpl(unsigned int                const  outputDimWrt,
                              unsigned int                const  inputDimWrt,
                              ref_vector<Eigen::VectorXd> const& input);

    virtual void ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                   unsigned int                const  inputDimWrt,
                                   ref_vector<Eigen::VectorXd> const& input,
                                   Eigen::VectorXd             const& vec);

    // virtual void ApplyHessianImpl(int                         const  outputDimWrt,
    //                       int                         const  inputDimWrt1,
    //                       int                         const  inputDimWrt2,
    //                       ref_vector<Eigen::VectorXd> const& input,
    //                       Eigen::VectorXd             const& sensitivity
    //                       Eigen::VectorXd             const& vec);


    // virtual void ApplyHessianByFD(int                         const  outputDimWrt,
    //                       int                         const  inputDimWrt1,
    //                       int                         const  inputDimWrt2,
    //                       ref_vector<Eigen::VectorXd> const& input,
    //                       Eigen::VectorXd             const& sensitivity
    //                       Eigen::VectorXd             const& vec);

  private:
    template<typename NextType, typename... Args>                                                                                               \
    inline Eigen::VectorXd const& GradientRecurse(unsigned int outWrt, unsigned int inWrt, ref_vector<Eigen::VectorXd>& vec, NextType const& ith, Args const&... args) {     \
        static_assert(std::is_same<Eigen::VectorXd, NextType>::value, "In ModPiece::Gradient, cannot cast input to Eigen::VectorXd."); \
        vec.push_back(std::cref((NextType&)ith));                                                                                               \
        return GradientRecurse(outWrt, inWrt, vec, args...);                                                                              \
    }

    template<typename NextType>                                                                                               \
    inline Eigen::VectorXd const& GradientRecurse(unsigned int outWrt, unsigned int inWrt, ref_vector<Eigen::VectorXd>& vec, NextType const& last, Eigen::VectorXd const& sens) {     \
        static_assert(std::is_same<Eigen::VectorXd, NextType>::value, "In ModPiece::Gradient, cannot cast input to Eigen::VectorXd."); \
        vec.push_back(std::cref((NextType&)last));                                                                                               \
        return Gradient(outWrt, inWrt, vec, sens);                                                                              \
    }

    template<typename... Args>                                                                                               \
    inline Eigen::MatrixXd const& Jacobian(unsigned int outWrt, unsigned int inWrt, ref_vector<Eigen::VectorXd>& vec, Eigen::VectorXd const& ith, Args const&... args) {     \
        vec.push_back(std::cref(ith));                                                                                               \
        return Jacobian(outWrt, inWrt, vec, args...);                                                                              \
    }
    inline Eigen::MatrixXd const& Jacobian(unsigned int outWrt, unsigned int inWrt, ref_vector<Eigen::VectorXd>& vec, Eigen::VectorXd const& last) {
      vec.push_back(std::cref(last));
      return Jacobian(outWrt, inWrt, vec);
    }

    template<typename... Args>                                                                                               \
    inline Eigen::MatrixXd JacobianByFD(unsigned int outWrt, unsigned int inWrt, ref_vector<Eigen::VectorXd>& vec, Eigen::VectorXd const& ith, Args const&... args) {     \
        vec.push_back(std::cref(ith));                                                                                               \
        return JacobianByFD(outWrt, inWrt, vec, args...);                                                                              \
    }
    inline Eigen::MatrixXd JacobianByFD(unsigned int outWrt, unsigned int inWrt, ref_vector<Eigen::VectorXd>& vec, Eigen::VectorXd const& last) {
      vec.push_back(std::cref(last));
      return JacobianByFD(outWrt, inWrt, vec);
    }

    template<typename NextType, typename... Args>
    inline Eigen::MatrixXd ApplyJacobianByFD(unsigned int outWrt, unsigned int  inWrt, ref_vector<Eigen::VectorXd>& vec, NextType const& ith, Args const&... args) {
      vec.push_back(std::cref(ith));
      return ApplyJacobianByFD(outWrt, inWrt, vec, args...);
    }
    template<typename NextType>
    inline Eigen::MatrixXd ApplyJacobianByFD(unsigned int outWrt, unsigned int inWrt, ref_vector<Eigen::VectorXd>& vec, NextType const& last, Eigen::VectorXd const& sens) {
      vec.push_back(std::cref(last));
      return ApplyJacobianByFD(outWrt, inWrt, vec, sens);
    }

  };

  }
}



#endif // #ifndef MODPIECE_H
