/*  File produced by Kranc */

#include "cctk.h"

extern "C" int TestSimpleWaveCaKernel_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "TestSimpleWaveCaKernel";
  CCTK_RegisterBanner(banner);
  return 0;
}
