/*  File produced by Kranc */

#include "cctk.h"

extern "C" int TestSimpleWave_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "TestSimpleWave";
  CCTK_RegisterBanner(banner);
  return 0;
}
