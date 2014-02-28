/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ConservationCalculation_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ConservationCalculation";
  CCTK_RegisterBanner(banner);
  return 0;
}
