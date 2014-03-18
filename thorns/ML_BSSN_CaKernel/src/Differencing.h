#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder21(u) ((-KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1o2dx)
#else
#  define PDstandardNthfdOrder21(u) (PDstandardNthfdOrder21_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1o2dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder22(u) ((-KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1o2dy)
#else
#  define PDstandardNthfdOrder22(u) (PDstandardNthfdOrder22_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1o2dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder23(u) ((-KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1))*p1o2dz)
#else
#  define PDstandardNthfdOrder23(u) (PDstandardNthfdOrder23_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder22_impl(u, p1o2dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder41(u) ((-8*KRANC_GFOFFSET3D(u,-1,0,0) + 8*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,-2,0,0) - KRANC_GFOFFSET3D(u,2,0,0))*p1o12dx)
#else
#  define PDstandardNthfdOrder41(u) (PDstandardNthfdOrder41_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-8*KRANC_GFOFFSET3D(u,-1,0,0) + 8*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,-2,0,0) - KRANC_GFOFFSET3D(u,2,0,0))*p1o12dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder42(u) ((-8*KRANC_GFOFFSET3D(u,0,-1,0) + 8*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,-2,0) - KRANC_GFOFFSET3D(u,0,2,0))*p1o12dy)
#else
#  define PDstandardNthfdOrder42(u) (PDstandardNthfdOrder42_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-8*KRANC_GFOFFSET3D(u,0,-1,0) + 8*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,-2,0) - KRANC_GFOFFSET3D(u,0,2,0))*p1o12dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder43(u) ((-8*KRANC_GFOFFSET3D(u,0,0,-1) + 8*KRANC_GFOFFSET3D(u,0,0,1) + KRANC_GFOFFSET3D(u,0,0,-2) - KRANC_GFOFFSET3D(u,0,0,2))*p1o12dz)
#else
#  define PDstandardNthfdOrder43(u) (PDstandardNthfdOrder43_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder42_impl(u, p1o12dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder211(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1odx2)
#else
#  define PDstandardNthfdOrder211(u) (PDstandardNthfdOrder211_impl(u,p1odx2,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder211_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder211_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0))*p1odx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder222(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1ody2)
#else
#  define PDstandardNthfdOrder222(u) (PDstandardNthfdOrder222_impl(u,p1ody2,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder222_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder222_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0))*p1ody2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder233(u) ((-2*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1))*p1odz2)
#else
#  define PDstandardNthfdOrder233(u) (PDstandardNthfdOrder233_impl(u,p1odz2,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder233_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder233_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder222_impl(u, p1odz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder411(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) + KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o12dx2)
#else
#  define PDstandardNthfdOrder411(u) (PDstandardNthfdOrder411_impl(u,pm1o12dx2,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder411_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder411_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) + KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o12dx2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder422(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) + KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o12dy2)
#else
#  define PDstandardNthfdOrder422(u) (PDstandardNthfdOrder422_impl(u,pm1o12dy2,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder422_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder422_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) + KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o12dy2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder433(u) ((30*KRANC_GFOFFSET3D(u,0,0,0) - 16*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) + KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2))*pm1o12dz2)
#else
#  define PDstandardNthfdOrder433(u) (PDstandardNthfdOrder433_impl(u,pm1o12dz2,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder433_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder433_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder422_impl(u, pm1o12dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder212(u) ((KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy)
#else
#  define PDstandardNthfdOrder212(u) (PDstandardNthfdOrder212_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder212_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder212_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder213(u) ((KRANC_GFOFFSET3D(u,-1,0,-1) - KRANC_GFOFFSET3D(u,-1,0,1) - KRANC_GFOFFSET3D(u,1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1))*p1o4dxdz)
#else
#  define PDstandardNthfdOrder213(u) (PDstandardNthfdOrder213_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder213_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder213_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder212_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder221(u) ((KRANC_GFOFFSET3D(u,-1,-1,0) - KRANC_GFOFFSET3D(u,-1,1,0) - KRANC_GFOFFSET3D(u,1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0))*p1o4dxdy)
#else
#  define PDstandardNthfdOrder221(u) (PDstandardNthfdOrder221_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder221_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder221_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder212_impl(u, p1o4dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder223(u) ((KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz)
#else
#  define PDstandardNthfdOrder223(u) (PDstandardNthfdOrder223_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder223_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder223_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder231(u) ((KRANC_GFOFFSET3D(u,-1,0,-1) - KRANC_GFOFFSET3D(u,-1,0,1) - KRANC_GFOFFSET3D(u,1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1))*p1o4dxdz)
#else
#  define PDstandardNthfdOrder231(u) (PDstandardNthfdOrder231_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder231_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder231_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder212_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder232(u) ((KRANC_GFOFFSET3D(u,0,-1,-1) - KRANC_GFOFFSET3D(u,0,-1,1) - KRANC_GFOFFSET3D(u,0,1,-1) + KRANC_GFOFFSET3D(u,0,1,1))*p1o4dydz)
#else
#  define PDstandardNthfdOrder232(u) (PDstandardNthfdOrder232_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder232_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder232_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder223_impl(u, p1o4dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder412(u) ((-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy)
#else
#  define PDstandardNthfdOrder412(u) (PDstandardNthfdOrder412_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder412_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder412_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder413(u) ((-64*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 64*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 8*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 8*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) + KRANC_GFOFFSET3D(u,-2,0,-2) - KRANC_GFOFFSET3D(u,-2,0,2) - KRANC_GFOFFSET3D(u,2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2))*p1o144dxdz)
#else
#  define PDstandardNthfdOrder413(u) (PDstandardNthfdOrder413_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder413_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder413_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder412_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder421(u) ((-64*(KRANC_GFOFFSET3D(u,-1,1,0) + KRANC_GFOFFSET3D(u,1,-1,0)) + 64*(KRANC_GFOFFSET3D(u,-1,-1,0) + KRANC_GFOFFSET3D(u,1,1,0)) + 8*(KRANC_GFOFFSET3D(u,-1,2,0) + KRANC_GFOFFSET3D(u,1,-2,0) + KRANC_GFOFFSET3D(u,-2,1,0) + KRANC_GFOFFSET3D(u,2,-1,0)) - 8*(KRANC_GFOFFSET3D(u,-1,-2,0) + KRANC_GFOFFSET3D(u,1,2,0) + KRANC_GFOFFSET3D(u,-2,-1,0) + KRANC_GFOFFSET3D(u,2,1,0)) + KRANC_GFOFFSET3D(u,-2,-2,0) - KRANC_GFOFFSET3D(u,-2,2,0) - KRANC_GFOFFSET3D(u,2,-2,0) + KRANC_GFOFFSET3D(u,2,2,0))*p1o144dxdy)
#else
#  define PDstandardNthfdOrder421(u) (PDstandardNthfdOrder421_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder421_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder421_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder412_impl(u, p1o144dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder423(u) ((-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz)
#else
#  define PDstandardNthfdOrder423(u) (PDstandardNthfdOrder423_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder423_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder423_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder431(u) ((-64*(KRANC_GFOFFSET3D(u,-1,0,1) + KRANC_GFOFFSET3D(u,1,0,-1)) + 64*(KRANC_GFOFFSET3D(u,-1,0,-1) + KRANC_GFOFFSET3D(u,1,0,1)) + 8*(KRANC_GFOFFSET3D(u,-1,0,2) + KRANC_GFOFFSET3D(u,1,0,-2) + KRANC_GFOFFSET3D(u,-2,0,1) + KRANC_GFOFFSET3D(u,2,0,-1)) - 8*(KRANC_GFOFFSET3D(u,-1,0,-2) + KRANC_GFOFFSET3D(u,1,0,2) + KRANC_GFOFFSET3D(u,-2,0,-1) + KRANC_GFOFFSET3D(u,2,0,1)) + KRANC_GFOFFSET3D(u,-2,0,-2) - KRANC_GFOFFSET3D(u,-2,0,2) - KRANC_GFOFFSET3D(u,2,0,-2) + KRANC_GFOFFSET3D(u,2,0,2))*p1o144dxdz)
#else
#  define PDstandardNthfdOrder431(u) (PDstandardNthfdOrder431_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder431_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder431_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder412_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder432(u) ((-64*(KRANC_GFOFFSET3D(u,0,-1,1) + KRANC_GFOFFSET3D(u,0,1,-1)) + 64*(KRANC_GFOFFSET3D(u,0,-1,-1) + KRANC_GFOFFSET3D(u,0,1,1)) + 8*(KRANC_GFOFFSET3D(u,0,-1,2) + KRANC_GFOFFSET3D(u,0,1,-2) + KRANC_GFOFFSET3D(u,0,-2,1) + KRANC_GFOFFSET3D(u,0,2,-1)) - 8*(KRANC_GFOFFSET3D(u,0,-1,-2) + KRANC_GFOFFSET3D(u,0,1,2) + KRANC_GFOFFSET3D(u,0,-2,-1) + KRANC_GFOFFSET3D(u,0,2,1)) + KRANC_GFOFFSET3D(u,0,-2,-2) - KRANC_GFOFFSET3D(u,0,-2,2) - KRANC_GFOFFSET3D(u,0,2,-2) + KRANC_GFOFFSET3D(u,0,2,2))*p1o144dydz)
#else
#  define PDstandardNthfdOrder432(u) (PDstandardNthfdOrder432_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL PDstandardNthfdOrder432_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDstandardNthfdOrder432_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder423_impl(u, p1o144dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder21(u) ((6*KRANC_GFOFFSET3D(u,0,0,0) - 4*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) + KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o16dx)
#else
#  define PDdissipationNthfdOrder21(u) (PDdissipationNthfdOrder21_impl(u,pm1o16dx,cdj,cdk))
static CCTK_REAL PDdissipationNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o16dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o16dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (6*KRANC_GFOFFSET3D(u,0,0,0) - 4*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) + KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o16dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder22(u) ((6*KRANC_GFOFFSET3D(u,0,0,0) - 4*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) + KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o16dy)
#else
#  define PDdissipationNthfdOrder22(u) (PDdissipationNthfdOrder22_impl(u,pm1o16dy,cdj,cdk))
static CCTK_REAL PDdissipationNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o16dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o16dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (6*KRANC_GFOFFSET3D(u,0,0,0) - 4*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) + KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o16dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder23(u) ((6*KRANC_GFOFFSET3D(u,0,0,0) - 4*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) + KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2))*pm1o16dz)
#else
#  define PDdissipationNthfdOrder23(u) (PDdissipationNthfdOrder23_impl(u,pm1o16dz,cdj,cdk))
static CCTK_REAL PDdissipationNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o16dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o16dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDdissipationNthfdOrder22_impl(u, pm1o16dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder41(u) ((-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 6*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o64dx)
#else
#  define PDdissipationNthfdOrder41(u) (PDdissipationNthfdOrder41_impl(u,p1o64dx,cdj,cdk))
static CCTK_REAL PDdissipationNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o64dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o64dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 6*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o64dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder42(u) ((-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 6*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o64dy)
#else
#  define PDdissipationNthfdOrder42(u) (PDdissipationNthfdOrder42_impl(u,p1o64dy,cdj,cdk))
static CCTK_REAL PDdissipationNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o64dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o64dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 6*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o64dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder43(u) ((-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) - 6*(KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2)) + KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3))*p1o64dz)
#else
#  define PDdissipationNthfdOrder43(u) (PDdissipationNthfdOrder43_impl(u,p1o64dz,cdj,cdk))
static CCTK_REAL PDdissipationNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o64dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDdissipationNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o64dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDdissipationNthfdOrder42_impl(u, p1o64dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder21(u) ((3*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,2*dir1,0,0) - 4*KRANC_GFOFFSET3D(u,dir1,0,0))*pm1o2dx*dir1)
#else
#  define PDupwindNthfdOrder21(u) (PDupwindNthfdOrder21_impl(u,pm1o2dx,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (3*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,2*dir1,0,0) - 4*KRANC_GFOFFSET3D(u,dir1,0,0))*pm1o2dx*dir1;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder41(u) ((-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,2*dir1,0,0) + KRANC_GFOFFSET3D(u,3*dir1,0,0) - 3*KRANC_GFOFFSET3D(u,-dir1,0,0) + 18*KRANC_GFOFFSET3D(u,dir1,0,0))*p1o12dx*dir1)
#else
#  define PDupwindNthfdOrder41(u) (PDupwindNthfdOrder41_impl(u,p1o12dx,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,2*dir1,0,0) + KRANC_GFOFFSET3D(u,3*dir1,0,0) - 3*KRANC_GFOFFSET3D(u,-dir1,0,0) + 18*KRANC_GFOFFSET3D(u,dir1,0,0))*p1o12dx*dir1;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder21(u) ((-4*KRANC_GFOFFSET3D(u,-1,0,0) + 4*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,-2,0,0) - KRANC_GFOFFSET3D(u,2,0,0))*p1o4dx)
#else
#  define PDupwindNthAntifdOrder21(u) (PDupwindNthAntifdOrder21_impl(u,p1o4dx,cdj,cdk))
static CCTK_REAL PDupwindNthAntifdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAntifdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-4*KRANC_GFOFFSET3D(u,-1,0,0) + 4*KRANC_GFOFFSET3D(u,1,0,0) + KRANC_GFOFFSET3D(u,-2,0,0) - KRANC_GFOFFSET3D(u,2,0,0))*p1o4dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder41(u) ((-21*KRANC_GFOFFSET3D(u,-1,0,0) + 21*KRANC_GFOFFSET3D(u,1,0,0) + 6*KRANC_GFOFFSET3D(u,-2,0,0) - 6*KRANC_GFOFFSET3D(u,2,0,0) - KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o24dx)
#else
#  define PDupwindNthAntifdOrder41(u) (PDupwindNthAntifdOrder41_impl(u,p1o24dx,cdj,cdk))
static CCTK_REAL PDupwindNthAntifdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o24dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAntifdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o24dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-21*KRANC_GFOFFSET3D(u,-1,0,0) + 21*KRANC_GFOFFSET3D(u,1,0,0) + 6*KRANC_GFOFFSET3D(u,-2,0,0) - 6*KRANC_GFOFFSET3D(u,2,0,0) - KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o24dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder21(u) ((6*KRANC_GFOFFSET3D(u,0,0,0) - 4*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) + KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o4dx)
#else
#  define PDupwindNthSymmfdOrder21(u) (PDupwindNthSymmfdOrder21_impl(u,pm1o4dx,cdj,cdk))
static CCTK_REAL PDupwindNthSymmfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o4dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymmfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o4dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (6*KRANC_GFOFFSET3D(u,0,0,0) - 4*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) + KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0))*pm1o4dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder41(u) ((-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 6*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o24dx)
#else
#  define PDupwindNthSymmfdOrder41(u) (PDupwindNthSymmfdOrder41_impl(u,p1o24dx,cdj,cdk))
static CCTK_REAL PDupwindNthSymmfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o24dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymmfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o24dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,-1,0,0) + KRANC_GFOFFSET3D(u,1,0,0)) - 6*(KRANC_GFOFFSET3D(u,-2,0,0) + KRANC_GFOFFSET3D(u,2,0,0)) + KRANC_GFOFFSET3D(u,-3,0,0) + KRANC_GFOFFSET3D(u,3,0,0))*p1o24dx;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided1(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,dir1,0,0))*p1odx*dir1)
#else
#  define PDonesided1(u) (PDonesided1_impl(u,p1odx,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odx, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,dir1,0,0))*p1odx*dir1;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder22(u) ((3*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,2*dir2,0) - 4*KRANC_GFOFFSET3D(u,0,dir2,0))*pm1o2dy*dir2)
#else
#  define PDupwindNthfdOrder22(u) (PDupwindNthfdOrder22_impl(u,pm1o2dy,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (3*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,2*dir2,0) - 4*KRANC_GFOFFSET3D(u,0,dir2,0))*pm1o2dy*dir2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder42(u) ((-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,0,2*dir2,0) + KRANC_GFOFFSET3D(u,0,3*dir2,0) - 3*KRANC_GFOFFSET3D(u,0,-dir2,0) + 18*KRANC_GFOFFSET3D(u,0,dir2,0))*p1o12dy*dir2)
#else
#  define PDupwindNthfdOrder42(u) (PDupwindNthfdOrder42_impl(u,p1o12dy,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,0,2*dir2,0) + KRANC_GFOFFSET3D(u,0,3*dir2,0) - 3*KRANC_GFOFFSET3D(u,0,-dir2,0) + 18*KRANC_GFOFFSET3D(u,0,dir2,0))*p1o12dy*dir2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder22(u) ((-4*KRANC_GFOFFSET3D(u,0,-1,0) + 4*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,-2,0) - KRANC_GFOFFSET3D(u,0,2,0))*p1o4dy)
#else
#  define PDupwindNthAntifdOrder22(u) (PDupwindNthAntifdOrder22_impl(u,p1o4dy,cdj,cdk))
static CCTK_REAL PDupwindNthAntifdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAntifdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-4*KRANC_GFOFFSET3D(u,0,-1,0) + 4*KRANC_GFOFFSET3D(u,0,1,0) + KRANC_GFOFFSET3D(u,0,-2,0) - KRANC_GFOFFSET3D(u,0,2,0))*p1o4dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder42(u) ((-21*KRANC_GFOFFSET3D(u,0,-1,0) + 21*KRANC_GFOFFSET3D(u,0,1,0) + 6*KRANC_GFOFFSET3D(u,0,-2,0) - 6*KRANC_GFOFFSET3D(u,0,2,0) - KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o24dy)
#else
#  define PDupwindNthAntifdOrder42(u) (PDupwindNthAntifdOrder42_impl(u,p1o24dy,cdj,cdk))
static CCTK_REAL PDupwindNthAntifdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o24dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAntifdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o24dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-21*KRANC_GFOFFSET3D(u,0,-1,0) + 21*KRANC_GFOFFSET3D(u,0,1,0) + 6*KRANC_GFOFFSET3D(u,0,-2,0) - 6*KRANC_GFOFFSET3D(u,0,2,0) - KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o24dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder22(u) ((6*KRANC_GFOFFSET3D(u,0,0,0) - 4*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) + KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o4dy)
#else
#  define PDupwindNthSymmfdOrder22(u) (PDupwindNthSymmfdOrder22_impl(u,pm1o4dy,cdj,cdk))
static CCTK_REAL PDupwindNthSymmfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o4dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymmfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o4dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (6*KRANC_GFOFFSET3D(u,0,0,0) - 4*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) + KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0))*pm1o4dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder42(u) ((-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 6*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o24dy)
#else
#  define PDupwindNthSymmfdOrder42(u) (PDupwindNthSymmfdOrder42_impl(u,p1o24dy,cdj,cdk))
static CCTK_REAL PDupwindNthSymmfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o24dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymmfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o24dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,0,-1,0) + KRANC_GFOFFSET3D(u,0,1,0)) - 6*(KRANC_GFOFFSET3D(u,0,-2,0) + KRANC_GFOFFSET3D(u,0,2,0)) + KRANC_GFOFFSET3D(u,0,-3,0) + KRANC_GFOFFSET3D(u,0,3,0))*p1o24dy;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided2(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,dir2,0))*p1ody*dir2)
#else
#  define PDonesided2(u) (PDonesided2_impl(u,p1ody,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1ody, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return (-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,dir2,0))*p1ody*dir2;
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder23(u) ((3*KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,0,2*dir3) - 4*KRANC_GFOFFSET3D(u,0,0,dir3))*pm1o2dz*dir3)
#else
#  define PDupwindNthfdOrder23(u) (PDupwindNthfdOrder23_impl(u,pm1o2dz,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthfdOrder22_impl(u, pm1o2dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder43(u) ((-10*KRANC_GFOFFSET3D(u,0,0,0) - 6*KRANC_GFOFFSET3D(u,0,0,2*dir3) + KRANC_GFOFFSET3D(u,0,0,3*dir3) - 3*KRANC_GFOFFSET3D(u,0,0,-dir3) + 18*KRANC_GFOFFSET3D(u,0,0,dir3))*p1o12dz*dir3)
#else
#  define PDupwindNthfdOrder43(u) (PDupwindNthfdOrder43_impl(u,p1o12dz,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDupwindNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthfdOrder42_impl(u, p1o12dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder23(u) ((-4*KRANC_GFOFFSET3D(u,0,0,-1) + 4*KRANC_GFOFFSET3D(u,0,0,1) + KRANC_GFOFFSET3D(u,0,0,-2) - KRANC_GFOFFSET3D(u,0,0,2))*p1o4dz)
#else
#  define PDupwindNthAntifdOrder23(u) (PDupwindNthAntifdOrder23_impl(u,p1o4dz,cdj,cdk))
static CCTK_REAL PDupwindNthAntifdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAntifdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o4dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthAntifdOrder22_impl(u, p1o4dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder43(u) ((-21*KRANC_GFOFFSET3D(u,0,0,-1) + 21*KRANC_GFOFFSET3D(u,0,0,1) + 6*KRANC_GFOFFSET3D(u,0,0,-2) - 6*KRANC_GFOFFSET3D(u,0,0,2) - KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3))*p1o24dz)
#else
#  define PDupwindNthAntifdOrder43(u) (PDupwindNthAntifdOrder43_impl(u,p1o24dz,cdj,cdk))
static CCTK_REAL PDupwindNthAntifdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o24dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthAntifdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o24dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthAntifdOrder42_impl(u, p1o24dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder23(u) ((6*KRANC_GFOFFSET3D(u,0,0,0) - 4*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) + KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2))*pm1o4dz)
#else
#  define PDupwindNthSymmfdOrder23(u) (PDupwindNthSymmfdOrder23_impl(u,pm1o4dz,cdj,cdk))
static CCTK_REAL PDupwindNthSymmfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o4dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymmfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL pm1o4dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthSymmfdOrder22_impl(u, pm1o4dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder43(u) ((-20*KRANC_GFOFFSET3D(u,0,0,0) + 15*(KRANC_GFOFFSET3D(u,0,0,-1) + KRANC_GFOFFSET3D(u,0,0,1)) - 6*(KRANC_GFOFFSET3D(u,0,0,-2) + KRANC_GFOFFSET3D(u,0,0,2)) + KRANC_GFOFFSET3D(u,0,0,-3) + KRANC_GFOFFSET3D(u,0,0,3))*p1o24dz)
#else
#  define PDupwindNthSymmfdOrder43(u) (PDupwindNthSymmfdOrder43_impl(u,p1o24dz,cdj,cdk))
static CCTK_REAL PDupwindNthSymmfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o24dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDupwindNthSymmfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1o24dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthSymmfdOrder42_impl(u, p1o24dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided3(u) ((-KRANC_GFOFFSET3D(u,0,0,0) + KRANC_GFOFFSET3D(u,0,0,dir3))*p1odz*dir3)
#else
#  define PDonesided3(u) (PDonesided3_impl(u,p1odz,cdj,cdk,dir1,dir2,dir3))
static CCTK_REAL PDonesided3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL PDonesided3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL p1odz, const ptrdiff_t cdj, const ptrdiff_t cdk, const ptrdiff_t dir1, const ptrdiff_t dir2, const ptrdiff_t dir3)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDonesided2_impl(u, p1odz, cdk, cdj);
}
#endif

