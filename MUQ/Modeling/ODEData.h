#ifndef ODEDATA_H_
#define ODEDATA_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Modeling {
    /// A helper class for Sundial's time integration
    /**
       Stores necessary muq::Modeling::WorkPiece's and thier inputs so Sundials can hold that information.
     */
    class ODEData {
    public:

      /// Construct basic ode data
      /**
	 @param[in] rhs A muq::Modeling::WorkPiece that evalautes the right hand side of an ODE
	 @param[in] inputs The inputs to the rhs --- the first is the state, the rest are constant in time
	 @param[in] wrtIn The input we are computing the derivative wrt --- negative indicates no derivative is being computed
	 @param[in] wrtOut The output we are computing the derivative of --- negative indicates no derivative is being computed
       */
      ODEData(std::shared_ptr<WorkPiece> rhs, ref_vector<boost::any> const& inputs, int const wrtIn, int const wrtOut);

      /// Construct with root function
      /**
	 @param[in] rhs A muq::Modeling::WorkPiece that evalautes the right hand side of an ODE
	 @param[in] root A muq::Modeling::WorkPiece that evalautes a function we are trying to find the root of along an orbit of the ODE
	 @param[in] inputs The inputs to the rhs --- the first is the state, the rest are constant in time
       */
      ODEData(std::shared_ptr<WorkPiece> rhs, std::shared_ptr<WorkPiece> root, ref_vector<boost::any> const& inputs);

      /// The right hand side of the ODE
      std::shared_ptr<WorkPiece> rhs;

      /// A function we are trying to find the root of along an orbit of the ODE (nullptr if we are not doing a root finding problem)
      std::shared_ptr<WorkPiece> root;

      /// The inputs to the rhs --- the first is the state, the rest are constant in time
      ref_vector<boost::any> inputs;

      /// The input we are computing the derivative wrt --- negative indicates no derivative is being computed
      int wrtIn = -1;

      /// The output we are computing the derivative of --- negative indicates no derivative is being computed
      int wrtOut = -1;
      
    private:
    };
  } // namespace Modeling
} // namespace muq

#endif
