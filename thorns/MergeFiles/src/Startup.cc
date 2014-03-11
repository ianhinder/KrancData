/*  File produced by Kranc */

#include "cctk.h"

extern "C" int MergeFiles_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "MergeFiles";
  CCTK_RegisterBanner(banner);
  return 0;
}
