/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ParamCheck_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ParamCheck";
  CCTK_RegisterBanner(banner);
  return 0;
}
