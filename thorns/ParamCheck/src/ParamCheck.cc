#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
extern "C" void ParamCheck_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (CCTK_MyProc(cctkGH) == 0)
  {
    if (aparam == 0)
    {
      CCTK_WARN(CCTK_WARN_ABORT, "aparam must be == 0");
    }
  }
}
