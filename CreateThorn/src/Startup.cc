/*  File produced by Kranc */

#include "cctk.h"

extern "C" int CreateThorn_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "CreateThorn";
  CCTK_RegisterBanner(banner);
  return 0;
}
