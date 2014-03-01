/*  File produced by Kranc */

#include "cctk.h"

extern "C" int LoopControlOne_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "LoopControlOne";
  CCTK_RegisterBanner(banner);
  return 0;
}
