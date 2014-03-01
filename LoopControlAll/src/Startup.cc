/*  File produced by Kranc */

#include "cctk.h"

extern "C" int LoopControlAll_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "LoopControlAll";
  CCTK_RegisterBanner(banner);
  return 0;
}
