/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_BSSN_NoVec_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ML_BSSN_NoVec";
  CCTK_RegisterBanner(banner);
  return 0;
}
