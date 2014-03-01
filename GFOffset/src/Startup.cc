/*  File produced by Kranc */

#include "cctk.h"

extern "C" int GFOffset_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "GFOffset";
  CCTK_RegisterBanner(banner);
  return 0;
}
