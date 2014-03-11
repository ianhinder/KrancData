/*  File produced by Kranc */

#include "cctk.h"

extern "C" int CountOperations_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "CountOperations";
  CCTK_RegisterBanner(banner);
  return 0;
}
