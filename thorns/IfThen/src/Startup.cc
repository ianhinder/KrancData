/*  File produced by Kranc */

#include "cctk.h"

extern "C" int IfThen_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "IfThen";
  CCTK_RegisterBanner(banner);
  return 0;
}
