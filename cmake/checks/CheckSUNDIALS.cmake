# make sure that the boost graph library is available
set(CMAKE_REQUIRED_LIBRARIES ${SUNDIALS_LIBRARIES} m)
set(CMAKE_REQUIRED_INCLUDES ${SUNDIALS_INCLUDE_DIRS})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-unresolved-symbols=ignore-in-shared-libs")

CHECK_CXX_SOURCE_COMPILES(
"
#include <idas/idas.h>           
#include <idas/idas_spgmr.h>     
#include <idas/idas_spbcgs.h>    /* prototypes & constants for CVSPBCG solver */
#include <idas/idas_sptfqmr.h>   /* prototypes & constants for SPTFQMR solver */
#include <idas/idas_dense.h>

#include <nvector/nvector_serial.h> 
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

static int idasJac(long int N,
                   realtype time,
				   realtype alpha,
                   N_Vector state,
				   N_Vector deriv,
                   N_Vector resid,
                   DlsMat   jac,
                   void    *user_data,
                   N_Vector tmp1,
                   N_Vector tmp2,
                   N_Vector tmp3){return 1;};
int main(){

void    *ida_mem;
IDADlsSetDenseJacFn(ida_mem, idasJac);

return 0;
}
"
SUNDIALS_IDA_COMPILES)

if(NOT SUNDIALS_IDA_COMPILES)
	set(SUNDIALS_TEST_FAIL 1)
else()
	set(SUNDIALS_TEST_FAIL 0)
endif()
