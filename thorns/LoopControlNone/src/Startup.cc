/*  File produced by Kranc */

#include "cctk.h"

extern "C" int LoopControlNone_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "LoopControlNone";
  CCTK_RegisterBanner(banner);
  return 0;
}
