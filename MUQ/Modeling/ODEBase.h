#ifndef ODEBASE_H_
#define ODEBASE_H_

#include "boost/property_tree/ptree.hpp"

#include "MUQ/Modeling/AnyAlgebra.h"
#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Modeling/ODEData.h"

namespace muq {
  namespace Modeling {

    /// A bass class to integrate ODE's
    class ODEBase : public WorkPiece {
    public:

      /**
	 The first input is the initial state (at \f$t=0\f$).  It is also the first input to the right hand side muq::Modeling::WorkPiece.  This must a N_Vector type (the vectors that Sundials uses).

	 The next set of inputs are the inputs to the right hand side muq::Modeling::WorkPiece.  If the right hand side input takes 2 inputs besides the state, these correspond to inputs 2 and 3 of the muq::Modeling::ODEBase.   Their types are known if the types are known by the rhs muq::Modeling::WorkPiece.

	 Additional inputs depend on which child inherits from muq::Modeling::ODEBase

	 Parameters options (for boost::property_tree::ptree):
	 <ol>
	 <li> The multistep method (<EM>ODESolver.MultistepMethod</EM>)
	 <ul>
	 <li> <EM>BDF</EM>: Backward differentiation formulas (BDF) linear mutlistep method (Default)
	 <li> <EM>Adams</EM>: Adams-Moulton linear multistep method (for nonstiff problems)
	 </ul>
	 <li> The linear solver (<EM>ODESolver.LinearSolver</EM>)
	 <ul>
	 <li> <EM>Dense</EM>: Dense linear system solve (default)
	 <li> <EM>SPGMR</EM>
	 <li> <EM>SPBCG</EM>
	 <li> <EM>SPTFQMR</EM>
	 </ul>
	 <li> The nonlinear solver (<EM>ODESolver.NonlinearSolver</EM>)
	 <ul>
	 <li> <EM>Newton</EM>: Nonlinear system solution through Newton's method (Default)
	 <li> <EM>Iter</EM>: Nonlinear system solution through functional iterations
	 </ul>
	 <li> The relative tolerance (<EM>ODESolver.RelativeTolerance</EM>)
	 <ul>
	 <li> Defaults to 1e-8
	 </ul>
	 <li> The absolute tolerance (<EM>ODESolver.AbsoluteTolerance</EM>)
	 <ul>
	 <li> Defaults to 1e-8
	 </ul>
	 <li> The maximum time step size (<EM>ODESolver.MaxStepSize</EM>)
	 <ul>
	 <li> Defaults to 1.0
	 </ul>
	 </ol>
	 @param[in] rhs The right hand side of the ODE
	 @param[in] pt A boost::property_tree::ptree with options/tolerances for the ODE integrator
	 @param[in] algebra A muq::Modeling::AnyAlgebra used to manipulate the state and input parameters
       */
      ODEBase(std::shared_ptr<WorkPiece> rhs, boost::property_tree::ptree const& pt, std::shared_ptr<AnyAlgebra> algebra);

      virtual ~ODEBase();
      
    protected:

      /// Are we computing the Jacobian, the action of the Jacobian, or the action of the Jacobian transpose
      enum DerivativeMode {
	/// The Jacobian
	Jac,
	/// The action of the Jacobian
	JacAction,
	/// The action of the Jacobian transpose
	JacTransAction
      };

      /// Check the return flag of a Sundials function
      /**
	 @param[in] flagvalue The value of the Sundials flag
	 @param[in] funcname The name of the Sundials function
	 @param[in] opt An option to determine how to check the flag, 0: check if flag is nullptr, 1: flag is an int, check if flag<0 (indicates Sundials error)
	 \return false: failure, true: success
       */
      bool CheckFlag(void* flagvalue, std::string const& funcname, unsigned int const opt) const;

      /// Initialize the state vector to the initial conditions
      /**
	 @param[out] state The state vector (to be initialized)
	 @param[in] ic The initial conditions 
	 @param[in] dim The dimension of the state
       */
      void InitializeState(N_Vector& state, boost::any const& ic, unsigned int const dim) const;

      /// Alloc memory and set up options for the Sundials solver
      /**
	 @param[in] cvode_mem The Sundials solver
	 @param[in] state The initial state
	 @param[in] data An object that holds the RHS inputs and can evaluate the RHS
       */
      void CreateSolverMemory(void* cvode_mem, N_Vector const& state, std::shared_ptr<ODEData> data) const;

      /// Copy the values of an N_Vector
      /**
	 @param[out] copy A new vector whose size and values are the same as orig
	 @param[in] orig An existing vector to be copied
       */
      void DeepCopy(N_Vector& copy, N_Vector const& orig) const;

      /// Initialize the derivative output (muq::Modeling::WorkPiece::jacobian, muq::Modeling::WorkPiece::jacobianAction, or muq::Modeling::WorkPiece::jacobianTransposeAction)
      /**
	 @param[in] ntimes The number of output times
	 @param[in] stateSize The size of the state
	 @param[in] paramSize The size of the input parameter we are differenating wrt
	 @param[in] mode Are we computing the Jacobian, Jacobian action, or Jacobian transpose action?	 
       */
      void InitializeDerivative(unsigned int const ntimes, unsigned int const stateSize, unsigned int const paramSize, DerivativeMode const& mode);

      /// Deal with Sundials errors
      /**
	 Sundials will call this function if it runs into a problem
	 @param[in] error_code Sundials error code
	 @param[in] module The name of the CVODES module reporting the error
	 @param[in] function The name of the function in which the error occured
	 @param[in] msg The error message
	 @param[in] user_data A pointer to an muq::Modeling::ODEData
       */
      static void ErrorHandler(int error_code, const char *module, const char *function, char *msg, void *user_data);

      /// Evaluate the right hand side
      /**
	 @param[in] time The current time 
	 @param[in] state The current state
	 @param[out] deriv The derivative of the state with respect to time
	 @param[in] user_data A pointer to an muq::Modeling::ODEData
       */
      static int EvaluateRHS(realtype time, N_Vector state, N_Vector deriv, void *user_data);

      /// Evaluate the Jacobian of the right hand side
      /**
	 @param[in] N
	 @param[in] time The current time 
	 @param[in] state The current state
	 @param[in] rhs The derivative of the state with respect to time
	 @param[out] jac The Jacobian of the right hand side with respect to the current state
	 @param[in] user_data A pointer to an muq::Modeling::ODEData
	 @param[in] tmp1
	 @param[in] tmp2
	 @param[in] tmp3
       */
      static int RHSJacobian(long int N, realtype time, N_Vector state, N_Vector rhs, DlsMat jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

      /// Evaluate the action of the Jacobian of the right hand side
      /**
	 @param[in] v The vector the Jacobian is acting on
	 @param[out] Jv The action of the Jacobian on v
	 @param[in] time The current time 
	 @param[in] state The current state
	 @param[in] rhs The derivative of the state with respect to time
	 @param[in] user_data A pointer to an muq::Modeling::ODEData
	 @param[in] tmp
       */
      static int RHSJacobianAction(N_Vector v, N_Vector Jv, realtype time, N_Vector state, N_Vector rhs, void *user_data, N_Vector tmp);

      /// Sundials uses this function to compute the derivative of the state at each timestep
      /**
	 @param[in] Ns The number of sensitivities
	 @param[in] time The current time
	 @param[in] y The current state
	 @param[in] ydot The derivative of the current state with respect to time 
	 @param[in] ys 
	 @param[in] ySdot The sensitivties 
	 @param[in] user_data A pointer to an muq::Modeling::ODEData
       */
      static int ForwardSensitivityRHS(int Ns, realtype time, N_Vector y, N_Vector ydot, N_Vector *ys, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2);

      /// Set up the solver for sensitivity information
      /**
	 @param[in] cvode_mem The Sundials solver 
	 @param[in] paramSize The size of the input parameter we are differenating wrt
	 @param[in,out] sensState This will become the 'current' Jacobian 
       */
      void SetUpSensitivity(void *cvode_mem, unsigned int const paramSize, N_Vector *sensState) const;

      /// The current index and size of each output vector
      /**
	 @params[in] outputTimes The times the user has asked for
	 \return A vector of pairs --- first: the current index of this output vector (starts at 0), second: the size of that output vector
       */
      std::vector<std::pair<unsigned int, unsigned int> > TimeIndices(ref_vector<boost::any> const& outputTimes);

      /// Compute the next time to integrate to
      /**
	 @param[out] nextTime first: the next time to integrate to, second: the output index
	 @param[out] timeIndices Each element corresponds to a vector of desired times, first: the current index of that vector, second: the size of that vector
	 @param[in] outputTimes We want the state at these times 
       */
      bool NextTime(std::pair<double, int>& nextTime, std::vector<std::pair<unsigned int, unsigned int> >& timeIndices, ref_vector<boost::any> const& outputTimes) const;

      /// Clear the outputs, jacobian, jacobianAction, and jacobianTransposeAction
      void ClearResults();

      /// Which linear solver should we use?
      enum LinearSolver {
	/// Dense solver
	Dense,
	/// SPGMR
	SPGMR,
	/// SPBCG
	SPBCG,
	/// SPTFQMR
	SPTFQMR
      };

      /// The right hand side of the ODE
      std::shared_ptr<WorkPiece> rhs;

      /// An algebra to manipulate the state and parameters
      std::shared_ptr<AnyAlgebra> algebra;

      /// Linear solver method
      const std::string linSolver;
      LinearSolver slvr;

      /// The relative tolerance
      const double reltol;

      /// The absolute tolerance
      const double abstol;

      /// The maximum time step size
      const double maxStepSize;

      /// Multistep method
      int multiStep;
      
      /// Nonlinear solver method
      int solveMethod;

      /// N_Vector type name
      const std::string N_VectorName = typeid(N_Vector).name();

      /// std::vector<N_Vector> type name
      const std::string stdvecN_VectorName = typeid(std::vector<N_Vector>).name();

      /// DlsMat type name
      const std::string DlsMatName = typeid(DlsMat).name();

      /// std::vector<DlsMat> type name
      const std::string stdvecDlsMatName = typeid(std::vector<DlsMat>).name();

    private:

      /// Set the input and output types based on the rhs muq::Modeling::WorkPiece's
      void SetInputOutputTypes();

      /// Clear the outputs
      void ClearOutputs();

      /// Clear the jacobian
      void ClearJacobian();

      /// Clear the jacobianAction
      void ClearJacobianAction();

      /// Clear the jacobianTransposeAction
      void ClearJacobianTransposeAction();
    };
  
  } // namespace Modeling
} // namespace muq

#endif
