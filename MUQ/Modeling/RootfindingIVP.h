#ifndef ROOTFINDINGIVP_H_
#define ROOTFINDINGIVP_H_

#include <cvodes/cvodes.h> // prototypes for CVODE fcts. and consts.

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Modeling {
    /// A rootfinding initial value problem --- find the root of a function along an orbit of an ODE
    class RootfindingIVP : public WorkPiece {
    public:

      /**
	 The first input is the initial state (at \f$t=0\f$).  It is also the first input to the both right hand side and the root muq::Modeling::WorkPiece.  The type must be the same in both sub-models (if it is known). 

	 The next set of inputs are the inputs to the right hand side muq::Modeling::WorkPiece.  If the right hand side input takes 2 inputs besides the state, these correspond to inputs 2 and 3 of the root finder.   Their types are known if the types are known by the rhs muq::Modeling::WorkPiece.

	 The next set of inputs are the inputs to the root muq::Modeling::WorkPiece.  If the root input takes 2 inputs besides the state, these correspond to inputs 4 and 5 of the root finder.  Their types are known if the types are known by the rhs muq::Modeling::WorkPiece.

	 An optional final input is a vector of times to save the state.  If the user gives muq::Modeling::RootfinderIVP this (optional) input it integrates until it finds a root, saving the state at each of these times.

	 The first output is the state at the root.  This output has the same type as the first input to the right hand side and the root (if either is known).  Note that the "the root" is the first state such that any one of the outputs from the root function is zero.  It is boost::none if no root was found.

	 The second output is the time at which the root (first output) was obtained.

	 The third output exists if the user gives muq::Modeling::RootfinderIVP times to save the state (optional final input).  This output is a vector --- std::vector<StateType> --- of states at the specified times.
	 @param[in] rhs A muq::Modeling::WorkPiece that evaluates the right hand side of the ODE
	 @param[in] root A muq::Modeling::WorkPiece whose outputs are double's --- we integrate the ODE until we find the first root of one of these outputs
       */
      RootfindingIVP(std::shared_ptr<WorkPiece> rhs, std::shared_ptr<WorkPiece> root);

      virtual ~RootfindingIVP();
      
    private:

      virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override;

      /// Set the input and output types based on the rhs and root muq::Modeling::WorkPiece's
      void SetInputOutputTypes();

      /// The right hand side of the ODE
      std::shared_ptr<WorkPiece> rhs;

      /// The root function
      std::shared_ptr<WorkPiece> root;
      
    };
  } // namespace Modeling
} // namespace muq

#endif
