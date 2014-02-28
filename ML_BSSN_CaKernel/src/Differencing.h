#include <assert.h>
#include "vectors.h"

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder21(u) (kmul(p1o2dx,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,-1,0,0))))
#else
#  define PDstandardNthfdOrder21(u) (PDstandardNthfdOrder21_impl(u,p1o2dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o2dx,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,-1,0,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder22(u) (kmul(p1o2dy,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,-1,0))))
#else
#  define PDstandardNthfdOrder22(u) (PDstandardNthfdOrder22_impl(u,p1o2dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o2dy,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,-1,0)));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder23(u) (kmul(p1o2dz,ksub(KRANC_GFOFFSET3D(u,0,0,1),KRANC_GFOFFSET3D(u,0,0,-1))))
#else
#  define PDstandardNthfdOrder23(u) (PDstandardNthfdOrder23_impl(u,p1o2dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder22_impl(u, p1o2dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder41(u) (kmul(p1o12dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-8),kmsub(KRANC_GFOFFSET3D(u,1,0,0),ToReal(8),KRANC_GFOFFSET3D(u,2,0,0))))))
#else
#  define PDstandardNthfdOrder41(u) (PDstandardNthfdOrder41_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o12dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-8),kmsub(KRANC_GFOFFSET3D(u,1,0,0),ToReal(8),KRANC_GFOFFSET3D(u,2,0,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder42(u) (kmul(p1o12dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-8),kmsub(KRANC_GFOFFSET3D(u,0,1,0),ToReal(8),KRANC_GFOFFSET3D(u,0,2,0))))))
#else
#  define PDstandardNthfdOrder42(u) (PDstandardNthfdOrder42_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o12dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-8),kmsub(KRANC_GFOFFSET3D(u,0,1,0),ToReal(8),KRANC_GFOFFSET3D(u,0,2,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder43(u) (kmul(p1o12dz,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-8),kmsub(KRANC_GFOFFSET3D(u,0,0,1),ToReal(8),KRANC_GFOFFSET3D(u,0,0,2))))))
#else
#  define PDstandardNthfdOrder43(u) (PDstandardNthfdOrder43_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder42_impl(u, p1o12dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder211(u) (kmul(p1odx2,kadd(KRANC_GFOFFSET3D(u,-1,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-2),KRANC_GFOFFSET3D(u,1,0,0)))))
#else
#  define PDstandardNthfdOrder211(u) (PDstandardNthfdOrder211_impl(u,p1odx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder211_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder211_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1odx2,kadd(KRANC_GFOFFSET3D(u,-1,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-2),KRANC_GFOFFSET3D(u,1,0,0))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder222(u) (kmul(p1ody2,kadd(KRANC_GFOFFSET3D(u,0,-1,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-2),KRANC_GFOFFSET3D(u,0,1,0)))))
#else
#  define PDstandardNthfdOrder222(u) (PDstandardNthfdOrder222_impl(u,p1ody2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder222_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder222_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1ody2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1ody2,kadd(KRANC_GFOFFSET3D(u,0,-1,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-2),KRANC_GFOFFSET3D(u,0,1,0))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder233(u) (kmul(p1odz2,kadd(KRANC_GFOFFSET3D(u,0,0,-1),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-2),KRANC_GFOFFSET3D(u,0,0,1)))))
#else
#  define PDstandardNthfdOrder233(u) (PDstandardNthfdOrder233_impl(u,p1odz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder233_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder233_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder222_impl(u, p1odz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder411(u) (kmul(pm1o12dx2,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-16),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(30)))))))
#else
#  define PDstandardNthfdOrder411(u) (PDstandardNthfdOrder411_impl(u,pm1o12dx2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder411_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder411_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o12dx2,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-16),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(30))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder422(u) (kmul(pm1o12dy2,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-16),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(30)))))))
#else
#  define PDstandardNthfdOrder422(u) (PDstandardNthfdOrder422_impl(u,pm1o12dy2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder422_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder422_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o12dy2,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-16),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(30))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder433(u) (kmul(pm1o12dz2,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kadd(KRANC_GFOFFSET3D(u,0,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(-16),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(30)))))))
#else
#  define PDstandardNthfdOrder433(u) (PDstandardNthfdOrder433_impl(u,pm1o12dz2,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder433_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder433_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder422_impl(u, pm1o12dz2, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder212(u) (kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(KRANC_GFOFFSET3D(u,1,1,0),kadd(KRANC_GFOFFSET3D(u,1,-1,0),KRANC_GFOFFSET3D(u,-1,1,0))))))
#else
#  define PDstandardNthfdOrder212(u) (PDstandardNthfdOrder212_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder212_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder212_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(KRANC_GFOFFSET3D(u,1,1,0),kadd(KRANC_GFOFFSET3D(u,1,-1,0),KRANC_GFOFFSET3D(u,-1,1,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder213(u) (kmul(p1o4dxdz,kadd(KRANC_GFOFFSET3D(u,-1,0,-1),ksub(KRANC_GFOFFSET3D(u,1,0,1),kadd(KRANC_GFOFFSET3D(u,1,0,-1),KRANC_GFOFFSET3D(u,-1,0,1))))))
#else
#  define PDstandardNthfdOrder213(u) (PDstandardNthfdOrder213_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder213_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder213_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder212_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder221(u) (kmul(p1o4dxdy,kadd(KRANC_GFOFFSET3D(u,-1,-1,0),ksub(KRANC_GFOFFSET3D(u,1,1,0),kadd(KRANC_GFOFFSET3D(u,1,-1,0),KRANC_GFOFFSET3D(u,-1,1,0))))))
#else
#  define PDstandardNthfdOrder221(u) (PDstandardNthfdOrder221_impl(u,p1o4dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder221_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder221_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder212_impl(u, p1o4dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder223(u) (kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(KRANC_GFOFFSET3D(u,0,1,1),kadd(KRANC_GFOFFSET3D(u,0,1,-1),KRANC_GFOFFSET3D(u,0,-1,1))))))
#else
#  define PDstandardNthfdOrder223(u) (PDstandardNthfdOrder223_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder223_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder223_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(KRANC_GFOFFSET3D(u,0,1,1),kadd(KRANC_GFOFFSET3D(u,0,1,-1),KRANC_GFOFFSET3D(u,0,-1,1)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder231(u) (kmul(p1o4dxdz,kadd(KRANC_GFOFFSET3D(u,-1,0,-1),ksub(KRANC_GFOFFSET3D(u,1,0,1),kadd(KRANC_GFOFFSET3D(u,1,0,-1),KRANC_GFOFFSET3D(u,-1,0,1))))))
#else
#  define PDstandardNthfdOrder231(u) (PDstandardNthfdOrder231_impl(u,p1o4dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder231_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder231_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder212_impl(u, p1o4dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder232(u) (kmul(p1o4dydz,kadd(KRANC_GFOFFSET3D(u,0,-1,-1),ksub(KRANC_GFOFFSET3D(u,0,1,1),kadd(KRANC_GFOFFSET3D(u,0,1,-1),KRANC_GFOFFSET3D(u,0,-1,1))))))
#else
#  define PDstandardNthfdOrder232(u) (PDstandardNthfdOrder232_impl(u,p1o4dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder232_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder232_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder223_impl(u, p1o4dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder412(u) (kmul(p1o144dxdy,kadd(KRANC_GFOFFSET3D(u,-2,-2,0),kadd(KRANC_GFOFFSET3D(u,2,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(64))),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0))))))))
#else
#  define PDstandardNthfdOrder412(u) (PDstandardNthfdOrder412_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder412_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder412_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o144dxdy,kadd(KRANC_GFOFFSET3D(u,-2,-2,0),kadd(KRANC_GFOFFSET3D(u,2,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(64))),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder413(u) (kmul(p1o144dxdz,kadd(KRANC_GFOFFSET3D(u,-2,0,-2),kadd(KRANC_GFOFFSET3D(u,2,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(64))),KRANC_GFOFFSET3D(u,2,0,-2)),KRANC_GFOFFSET3D(u,-2,0,2))))))))
#else
#  define PDstandardNthfdOrder413(u) (PDstandardNthfdOrder413_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder413_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder413_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder412_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder421(u) (kmul(p1o144dxdy,kadd(KRANC_GFOFFSET3D(u,-2,-2,0),kadd(KRANC_GFOFFSET3D(u,2,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),ToReal(64))),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0))))))))
#else
#  define PDstandardNthfdOrder421(u) (PDstandardNthfdOrder421_impl(u,p1o144dxdy,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder421_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder421_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder412_impl(u, p1o144dxdy, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder423(u) (kmul(p1o144dydz,kadd(KRANC_GFOFFSET3D(u,0,-2,-2),kadd(KRANC_GFOFFSET3D(u,0,2,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(64))),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2))))))))
#else
#  define PDstandardNthfdOrder423(u) (PDstandardNthfdOrder423_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder423_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder423_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o144dydz,kadd(KRANC_GFOFFSET3D(u,0,-2,-2),kadd(KRANC_GFOFFSET3D(u,0,2,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(64))),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder431(u) (kmul(p1o144dxdz,kadd(KRANC_GFOFFSET3D(u,-2,0,-2),kadd(KRANC_GFOFFSET3D(u,2,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),ToReal(64))),KRANC_GFOFFSET3D(u,2,0,-2)),KRANC_GFOFFSET3D(u,-2,0,2))))))))
#else
#  define PDstandardNthfdOrder431(u) (PDstandardNthfdOrder431_impl(u,p1o144dxdz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder431_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder431_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder412_impl(u, p1o144dxdz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDstandardNthfdOrder432(u) (kmul(p1o144dydz,kadd(KRANC_GFOFFSET3D(u,0,-2,-2),kadd(KRANC_GFOFFSET3D(u,0,2,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),ToReal(-64),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),ToReal(-8),ksub(ksub(kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),ToReal(8),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),ToReal(64))),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2))))))))
#else
#  define PDstandardNthfdOrder432(u) (PDstandardNthfdOrder432_impl(u,p1o144dydz,cdj,cdk))
static CCTK_REAL_VEC PDstandardNthfdOrder432_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDstandardNthfdOrder432_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDstandardNthfdOrder423_impl(u, p1o144dydz, cdj, cdk);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder21(u) (kmul(pm1o16dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6)))))))
#else
#  define PDdissipationNthfdOrder21(u) (PDdissipationNthfdOrder21_impl(u,pm1o16dx,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o16dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o16dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o16dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder22(u) (kmul(pm1o16dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6)))))))
#else
#  define PDdissipationNthfdOrder22(u) (PDdissipationNthfdOrder22_impl(u,pm1o16dy,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o16dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o16dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o16dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder23(u) (kmul(pm1o16dz,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kadd(KRANC_GFOFFSET3D(u,0,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6)))))))
#else
#  define PDdissipationNthfdOrder23(u) (PDdissipationNthfdOrder23_impl(u,pm1o16dz,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o16dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o16dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDdissipationNthfdOrder22_impl(u, pm1o16dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder41(u) (kmul(p1o64dx,kadd(KRANC_GFOFFSET3D(u,-3,0,0),kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(15))))))))
#else
#  define PDdissipationNthfdOrder41(u) (PDdissipationNthfdOrder41_impl(u,p1o64dx,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o64dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o64dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o64dx,kadd(KRANC_GFOFFSET3D(u,-3,0,0),kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder42(u) (kmul(p1o64dy,kadd(KRANC_GFOFFSET3D(u,0,-3,0),kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(15))))))))
#else
#  define PDdissipationNthfdOrder42(u) (PDdissipationNthfdOrder42_impl(u,p1o64dy,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o64dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o64dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o64dy,kadd(KRANC_GFOFFSET3D(u,0,-3,0),kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDdissipationNthfdOrder43(u) (kmul(p1o64dz,kadd(KRANC_GFOFFSET3D(u,0,0,-3),kadd(KRANC_GFOFFSET3D(u,0,0,3),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(15))))))))
#else
#  define PDdissipationNthfdOrder43(u) (PDdissipationNthfdOrder43_impl(u,p1o64dz,cdj,cdk))
static CCTK_REAL_VEC PDdissipationNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o64dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDdissipationNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o64dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDdissipationNthfdOrder42_impl(u, p1o64dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder21(u) (kmul(pm1o2dx,kmul(dir1,kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(KRANC_GFOFFSET3D(u,1,0,0),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3)))))))
#else
#  define PDupwindNthfdOrder21(u) (PDupwindNthfdOrder21_impl(u,pm1o2dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder41(u) (kmul(p1o12dx,kmul(dir1,kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-10),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-6),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-3),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(18)))))))))
#else
#  define PDupwindNthfdOrder41(u) (PDupwindNthfdOrder41_impl(u,p1o12dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder21(u) (kmul(p1o4dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-4),kmsub(KRANC_GFOFFSET3D(u,1,0,0),ToReal(4),KRANC_GFOFFSET3D(u,2,0,0))))))
#else
#  define PDupwindNthAntifdOrder21(u) (PDupwindNthAntifdOrder21_impl(u,p1o4dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o4dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-4),kmsub(KRANC_GFOFFSET3D(u,1,0,0),ToReal(4),KRANC_GFOFFSET3D(u,2,0,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder41(u) (kmul(p1o24dx,kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-21),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-6),ksub(kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(6),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(21))),KRANC_GFOFFSET3D(u,-3,0,0)))))))
#else
#  define PDupwindNthAntifdOrder41(u) (PDupwindNthAntifdOrder41_impl(u,p1o24dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o24dx,kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,-1,0,0),ToReal(-21),kmadd(KRANC_GFOFFSET3D(u,2,0,0),ToReal(-6),ksub(kmadd(KRANC_GFOFFSET3D(u,-2,0,0),ToReal(6),kmul(KRANC_GFOFFSET3D(u,1,0,0),ToReal(21))),KRANC_GFOFFSET3D(u,-3,0,0))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder21(u) (kmul(pm1o4dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6)))))))
#else
#  define PDupwindNthSymmfdOrder21(u) (PDupwindNthSymmfdOrder21_impl(u,pm1o4dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o4dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o4dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o4dx,kadd(KRANC_GFOFFSET3D(u,-2,0,0),kadd(KRANC_GFOFFSET3D(u,2,0,0),kmadd(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder41(u) (kmul(p1o24dx,kadd(KRANC_GFOFFSET3D(u,-3,0,0),kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(15))))))))
#else
#  define PDupwindNthSymmfdOrder41(u) (PDupwindNthSymmfdOrder41_impl(u,p1o24dx,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder41_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o24dx,kadd(KRANC_GFOFFSET3D(u,-3,0,0),kadd(KRANC_GFOFFSET3D(u,3,0,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided1(u) (kmul(p1odx,kmul(dir1,ksub(KRANC_GFOFFSET3D(u,1,0,0),KRANC_GFOFFSET3D(u,0,0,0)))))
#else
#  define PDonesided1(u) (PDonesided1_impl(u,p1odx,cdj,cdk))
static CCTK_REAL_VEC PDonesided1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesided1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odx, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder22(u) (kmul(pm1o2dy,kmul(dir2,kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(KRANC_GFOFFSET3D(u,0,1,0),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3)))))))
#else
#  define PDupwindNthfdOrder22(u) (PDupwindNthfdOrder22_impl(u,pm1o2dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder42(u) (kmul(p1o12dy,kmul(dir2,kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-10),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-6),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-3),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(18)))))))))
#else
#  define PDupwindNthfdOrder42(u) (PDupwindNthfdOrder42_impl(u,p1o12dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder22(u) (kmul(p1o4dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-4),kmsub(KRANC_GFOFFSET3D(u,0,1,0),ToReal(4),KRANC_GFOFFSET3D(u,0,2,0))))))
#else
#  define PDupwindNthAntifdOrder22(u) (PDupwindNthAntifdOrder22_impl(u,p1o4dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o4dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-4),kmsub(KRANC_GFOFFSET3D(u,0,1,0),ToReal(4),KRANC_GFOFFSET3D(u,0,2,0)))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder42(u) (kmul(p1o24dy,kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-21),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-6),ksub(kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(6),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(21))),KRANC_GFOFFSET3D(u,0,-3,0)))))))
#else
#  define PDupwindNthAntifdOrder42(u) (PDupwindNthAntifdOrder42_impl(u,p1o24dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o24dy,kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,-1,0),ToReal(-21),kmadd(KRANC_GFOFFSET3D(u,0,2,0),ToReal(-6),ksub(kmadd(KRANC_GFOFFSET3D(u,0,-2,0),ToReal(6),kmul(KRANC_GFOFFSET3D(u,0,1,0),ToReal(21))),KRANC_GFOFFSET3D(u,0,-3,0))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder22(u) (kmul(pm1o4dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6)))))))
#else
#  define PDupwindNthSymmfdOrder22(u) (PDupwindNthSymmfdOrder22_impl(u,pm1o4dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o4dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o4dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(pm1o4dy,kadd(KRANC_GFOFFSET3D(u,0,-2,0),kadd(KRANC_GFOFFSET3D(u,0,2,0),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder42(u) (kmul(p1o24dy,kadd(KRANC_GFOFFSET3D(u,0,-3,0),kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(15))))))))
#else
#  define PDupwindNthSymmfdOrder42(u) (PDupwindNthSymmfdOrder42_impl(u,p1o24dy,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder42_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dy, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return kmul(p1o24dy,kadd(KRANC_GFOFFSET3D(u,0,-3,0),kadd(KRANC_GFOFFSET3D(u,0,3,0),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),ToReal(15)))))));
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided2(u) (kmul(p1ody,kmul(dir2,ksub(KRANC_GFOFFSET3D(u,0,1,0),KRANC_GFOFFSET3D(u,0,0,0)))))
#else
#  define PDonesided2(u) (PDonesided2_impl(u,p1ody,cdj,cdk))
static CCTK_REAL_VEC PDonesided2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1ody, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesided2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1ody, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder23(u) (kmul(pm1o2dz,kmul(dir3,kadd(KRANC_GFOFFSET3D(u,0,0,2),kmadd(KRANC_GFOFFSET3D(u,0,0,1),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(3)))))))
#else
#  define PDupwindNthfdOrder23(u) (PDupwindNthfdOrder23_impl(u,pm1o2dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o2dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthfdOrder43(u) (kmul(p1o12dz,kmul(dir3,kadd(KRANC_GFOFFSET3D(u,0,0,3),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-10),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-6),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-3),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(18)))))))))
#else
#  define PDupwindNthfdOrder43(u) (PDupwindNthfdOrder43_impl(u,p1o12dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder23(u) (kmul(p1o4dz,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-4),kmsub(KRANC_GFOFFSET3D(u,0,0,1),ToReal(4),KRANC_GFOFFSET3D(u,0,0,2))))))
#else
#  define PDupwindNthAntifdOrder23(u) (PDupwindNthAntifdOrder23_impl(u,p1o4dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o4dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthAntifdOrder22_impl(u, p1o4dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthAntifdOrder43(u) (kmul(p1o24dz,kadd(KRANC_GFOFFSET3D(u,0,0,3),kmadd(KRANC_GFOFFSET3D(u,0,0,-1),ToReal(-21),kmadd(KRANC_GFOFFSET3D(u,0,0,2),ToReal(-6),ksub(kmadd(KRANC_GFOFFSET3D(u,0,0,-2),ToReal(6),kmul(KRANC_GFOFFSET3D(u,0,0,1),ToReal(21))),KRANC_GFOFFSET3D(u,0,0,-3)))))))
#else
#  define PDupwindNthAntifdOrder43(u) (PDupwindNthAntifdOrder43_impl(u,p1o24dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthAntifdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthAntifdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthAntifdOrder42_impl(u, p1o24dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder23(u) (kmul(pm1o4dz,kadd(KRANC_GFOFFSET3D(u,0,0,-2),kadd(KRANC_GFOFFSET3D(u,0,0,2),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(-4),kmul(KRANC_GFOFFSET3D(u,0,0,0),ToReal(6)))))))
#else
#  define PDupwindNthSymmfdOrder23(u) (PDupwindNthSymmfdOrder23_impl(u,pm1o4dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o4dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o4dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthSymmfdOrder22_impl(u, pm1o4dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDupwindNthSymmfdOrder43(u) (kmul(p1o24dz,kadd(KRANC_GFOFFSET3D(u,0,0,-3),kadd(KRANC_GFOFFSET3D(u,0,0,3),kmadd(KRANC_GFOFFSET3D(u,0,0,0),ToReal(-20),kmadd(kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2)),ToReal(-6),kmul(kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),ToReal(15))))))))
#else
#  define PDupwindNthSymmfdOrder43(u) (PDupwindNthSymmfdOrder43_impl(u,p1o24dz,cdj,cdk))
static CCTK_REAL_VEC PDupwindNthSymmfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDupwindNthSymmfdOrder43_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o24dz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);
  return PDupwindNthSymmfdOrder42_impl(u, p1o24dz, cdk, cdj);
}
#endif

#ifndef KRANC_DIFF_FUNCTIONS
#  define PDonesided3(u) (kmul(p1odz,kmul(dir3,ksub(KRANC_GFOFFSET3D(u,0,0,1),KRANC_GFOFFSET3D(u,0,0,0)))))
#else
#  define PDonesided3(u) (PDonesided3_impl(u,p1odz,cdj,cdk))
static CCTK_REAL_VEC PDonesided3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;
static CCTK_REAL_VEC PDonesided3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1odz, const ptrdiff_t cdj, const ptrdiff_t cdk)
{ assert(0); return ToReal(1e30); /* ERROR */ }
#endif

