/*  File produced by Kranc */

#include "cctk.h"

extern "C" int TestSimpleWaveVectors_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "TestSimpleWaveVectors";
  CCTK_RegisterBanner(banner);
  return 0;
}
