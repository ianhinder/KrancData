/*  File produced by Kranc */

#include "cctk.h"

extern "C" int Analysis_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "Analysis";
  CCTK_RegisterBanner(banner);
  return 0;
}
