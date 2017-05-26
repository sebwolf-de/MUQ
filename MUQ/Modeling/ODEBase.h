#ifndef ODEBASE_H_
#define ODEBASE_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Modeling {

    class ODEBase : public WorkPiece {
    public:
      
      ODEBase(std::shared_ptr<WorkPiece> rhs);

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

      /// The right hand side of the ODE
      std::shared_ptr<WorkPiece> rhs;
      
    private:

      /// Set the input and output types based on the rhs muq::Modeling::WorkPiece's
      void SetInputOutputTypes();
    };
  
  } // namespace Modeling
} // namespace muq

#endif
