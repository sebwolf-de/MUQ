#ifndef ODE_H_
#define ODE_H_

#include "MUQ/Modeling/ODEBase.h"

namespace muq {
  namespace Modeling {
    class ODE : public ODEBase {
    public:
      /**
	 The first input is the initial state (at \f$t=0\f$).  It is also the first input to the right hand side muq::Modeling::WorkPiece.  This must a N_Vector type (the vectors that Sundials uses).

	 The next set of inputs are the inputs to the right hand side muq::Modeling::WorkPiece.  If the right hand side input takes 2 inputs besides the state, these correspond to inputs 2 and 3 of the muq::Modeling::ODEBase.   Their types are known if the types are known by the rhs muq::Modeling::WorkPiece.

	 Any inputs after the right hand sides inputs are either doubles or vectors of doubles.  The output type is the state at these times.  For example, if the last three inputs (after the RHS inputs) were 2.0 and std::vector<double>({1.0, 2.0}) and 1.5 then there would be three outputs --- the state at time 2.0, a vector of the state at times 1.0 and 2.0, and the state at time 1.5.
	 @param[in] rhs The right hand side of the ODE
	 @param[in] pt A boost::property_tree::ptree with options/tolerances for the ODE integrator
       */
      ODE(std::shared_ptr<ModPiece> rhs,  boost::property_tree::ptree const& pt);

      virtual ~ODE();

    private:

      /// Integrate the ODE forward in time
      /**
	 \f$M\f$ inputs:
	 <ul>
	 <li> The first \f$N\f$ are the inputs to muq::Modeling::ODE::rhs
	 <li> The second \fM-N\f$ are the times where we want to return the state (either vectors or scalars)
	 </ul>
	 @param[in] inputs The inputs (see description)
       */
      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      /// Evaluate the Jacobian of the state wrt each parameter
      /**
	 Returns the Jacobian of the state with respect to an input parameter at each timestep.  Returns a vector of Jacobians at each timestep if the output has more than one time.

	 \f$M\f$ inputs:
	 <ul>
	 <li> The first \f$N\f$ are the inputs to muq::Modeling::ODE::rhs
	 <li> The second \fM-N\f$ are the times where we want to return the state (either vectors or scalars)
	 </ul>
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
	 @param[in] wrtOut We are computing the derivative of this output
	 @param[in] inputs The inputs (see description)
       */
      virtual void JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<Eigen::VectorXd> const& inputs) override;

      /// Evaluate the action of the Jacobian on a given vector wrt each parameter
      /**
	 Returns the action of the Jacobian on a given vector with respect to an input parameter at each timestep.  Returns a vector of Jacobians actions at each timestep if the output has more than one time.

	 \f$M\f$ inputs:
	 <ul>
	 <li> The first \f$N\f$ are the inputs to muq::Modeling::ODE::rhs
	 <li> The second \fM-N\f$ are the times where we want to return the state (either vectors or scalars)
	 </ul>
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
	 @param[in] wrtOut We are computing the derivative of this output
   @param[in] inputs The inputs (see description)
	 @param[in] vec The vector the Jacobian is acting on at eact timestep
       */
      virtual void ApplyJacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& vec) override;

      /// Evaluate the action of the Jacobian transpose on a given vector wrt each parameter
      /**
	 Returns the action of the Jacobian transpose on a given vector with respect to an input parameter at each timestep.  Returns a vector of Jacobian transpose actions at each timestep if the output has more than one time.

	 \f$M\f$ inputs:
	 <ul>
	 <li> The first \f$N\f$ are the inputs to muq::Modeling::ODE::rhs
	 <li> The second \fM-N\f$ are the times where we want to return the state (either vectors or scalars)
	 </ul>
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
	 @param[in] wrtOut We are computing the derivative of this output
   @param[in] inputs The inputs (see description)
	 @param[in] vec The vector the Jacobian transpose is acting on at eact timestep
       */
      virtual void GradientImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<Eigen::VectorXd> const& inputs, Eigen::VectorXd const& vec) override;

      /// Integrate forward in time, possibly keeping track fo forward sensitivities
      /**
	 @param[in] inputs The inputs to the right hand side and the output times
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
	 @param[in] wrtOut We are computing the derivative of this output
	 @param[in] vec The vector the Jacobian transpose is acting on at eact timestep
	 @param[in] mode Are we computing the Jacobian, Jacobian action, or Jacobian transpose action?
       */
      void Integrate(ref_vector<Eigen::VectorXd> const& inputs, int const wrtIn = -1, int const wrtOut = -1,  N_Vector const& vec = nullptr, DerivativeMode const& mode = DerivativeMode::Jac);

      /// Integrate forward in time without computing derivative information
      /**
	 @param[in] cvode_mem The Sundials solver
	 @param[in,out] state The state --- begins at initial state, ends at final state
	 @param[in] outputTimes The output times the user has requested
       */
      void ForwardIntegration(void *cvode_mem, N_Vector& state, ref_vector<Eigen::VectorXd> const& outputTimes);

      /// Integrate forward in time, computing derivative information
      /**
	 @param[in] cvode_mem The Sundials solver
	 @param[in,out] state The state --- begins at initial state, ends at final state
	 @param[in] paramSize The size of the input parameter we are differenating wrt
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
	 @param[in] outputTimes The output times corresponding to the output we are differentiating
	 @param[in] rhsInputs The inputs to the right hand side
	 @param[in] mode Are we computing the Jacobian, Jacobian action, or Jacobian transpose action?
	 @param[in] vec The vector the Jacobian or its transpose is acting on at eact timestep
       */
      void ForwardSensitivity(void *cvode_mem, N_Vector& state, unsigned int const paramSize, unsigned int const wrtIn, boost::any const& outputTimes, ref_vector<Eigen::VectorXd> const& rhsInputs, DerivativeMode const& mode, N_Vector const& vec);

      /// Save the derivative information
      /**
	 @param[in] ntimes The number of output times
	 @param[in] timeIndex Index of the output vector to save the current derivative information in
	 @param[in] paramSize The size of the input parameter we are differenating wrt
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
	 @param[in] sensState The Jacobian of the state (resulting from Sundials solver) --- nullptr if we are not differentiating with respect to an input to the right hand side
	 @param[in] state The state --- begins at initial state, ends at final state
	 @param[in] rhsInputs The inputs to the right hand side
	 @param[in] vec The vector the Jacobian or its transpose is acting on at eact timestep
	 @param[in] mode Are we computing the Jacobian, Jacobian action, or Jacobian transpose action?
       */
      void SaveDerivative(unsigned int const ntimes, unsigned int const timeIndex, unsigned int const paramSize, unsigned int const wrtIn, N_Vector* sensState, N_Vector const& state, ref_vector<boost::any> rhsInputs, N_Vector const& vec, DerivativeMode const& mode);

      /// Save the Jacobian
      /**
	 @param[out] jac The Jacobian at this time step
	 @param[in] ncols The number of columns in the Jacobian
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
	 @param[in] sensState The Jacobian of the state (resulting from Sundials solver) --- nullptr if we are not differentiating with respect to an input to the right hand side
	 @param[in] state The state --- begins at initial state, ends at final state
	 @param[in] rhsInputs The inputs to the right hand side
       */
      void SaveJacobian(DlsMat& jac, unsigned int const ncols, unsigned int const wrtIn, N_Vector* sensState, N_Vector const& state, ref_vector<boost::any> rhsInputs) const;

      /// Save the action of the Jacobian
      /**
	 @param[out] jacAct The action of the Jacobian at this time step
	 @param[in] ncols The number of columns in the Jacobian
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
	 @param[in] sensState The Jacobian of the state (resulting from Sundials solver) --- nullptr if we are not differentiating with respect to an input to the right hand side
	 @param[in] state The state --- begins at initial state, ends at final state
	 @param[in] rhsInputs The inputs to the right hand side
	 @param[in] vec The vector the Jacobian is acting on at eact timestep
       */
      void SaveJacobianAction(N_Vector& jacAct, unsigned int const ncols, unsigned int const wrtIn, N_Vector* sensState, N_Vector const& state, ref_vector<boost::any> rhsInputs, N_Vector const& vec) const;

      /// Save the action of the Jacobian transpose
      /**
	 @param[out] jacTransAct The action of the Jacobian transpose at this time step
	 @param[in] ncols The number of columns in the Jacobian
	 @param[in] wrtIn We are computing the derivative with repsect to this parameter
	 @param[in] sensState The Jacobian of the state (resulting from Sundials solver) --- nullptr if we are not differentiating with respect to an input to the right hand side
	 @param[in] state The state --- begins at initial state, ends at final state
	 @param[in] rhsInputs The inputs to the right hand side
	 @param[in] vec The vector the Jacobian transpose is acting on at eact timestep
       */
      void SaveJacobianTransposeAction(N_Vector& jacTransAct, unsigned int const ncols, unsigned int const wrtIn, N_Vector* sensState, N_Vector const& state, ref_vector<boost::any> rhsInputs, N_Vector const& vec) const;
    };
  } // namespace Modeling
} // namespace muq

#endif
