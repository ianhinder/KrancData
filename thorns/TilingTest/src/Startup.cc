/*  File produced by Kranc */

#include "cctk.h"

extern "C" int TilingTest_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "TilingTest";
  CCTK_RegisterBanner(banner);
  return 0;
}
