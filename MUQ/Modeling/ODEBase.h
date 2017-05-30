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

      void CreateSolverMemory(void* cvode_mem, N_Vector const& state, std::shared_ptr<ODEData> data) const;

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

      /// Deal with Sundials errors
      /**
	 Sundials will call this function if it runs into a problem
	 @param[in] error_code Sundials error code
	 @param[in] module The name of the CVODES module reporting the error
	 @param[in] function The name of the function in which the error occured
	 @param[in] msg The error message
	 @param[in] eh_data A pointer to user data
       */
      static void ErrorHandler(int error_code, const char *module, const char *function, char *msg, void *eh_data);

      static int EvaluateRHS(realtype time, N_Vector state, N_Vector deriv, void *user_data);

      static int RHSJacobian(long int N, realtype time, N_Vector state, N_Vector rhs, DlsMat jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

      static int RHSJacobianAction(N_Vector v, N_Vector Jv, realtype time, N_Vector state, N_Vector rhs, void *user_data, N_Vector tmp);
      
    private:

      /// Set the input and output types based on the rhs muq::Modeling::WorkPiece's
      void SetInputOutputTypes();
    };
  
  } // namespace Modeling
} // namespace muq

#endif
