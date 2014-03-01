/*  File produced by Kranc */

#include "cctk.h"

extern "C" int TestSimpleWaveODE_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "TestSimpleWaveODE";
  CCTK_RegisterBanner(banner);
  return 0;
}
