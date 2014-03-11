/*  File produced by Kranc */

#define KRANC_C

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Kranc.hh"
#include "Differencing.h"
#include "loopcontrol.h"
#include "vectors.h"

/* Define macros used in calculations */
#define INITVALUE (42)
#define ScalarINV(x) ((CCTK_REAL)1.0 / (x))
#define ScalarSQR(x) ((x) * (x))
#define ScalarCUB(x) ((x) * ScalarSQR(x))
#define ScalarQAD(x) (ScalarSQR(ScalarSQR(x)))
#define INV(x) (kdiv(ToReal(1.0),x))
#define SQR(x) (kmul(x,x))
#define CUB(x) (kmul(x,SQR(x)))
#define QAD(x) (SQR(SQR(x)))

namespace ML_BSSN_DGFE {

extern "C" void ML_BSSN_DGFE_RHSStaticBoundary_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_DGFE_RHSStaticBoundary_calc_every != ML_BSSN_DGFE_RHSStaticBoundary_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_curvrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_dtlapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_dtlapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_dtshiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_Gammarhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_log_confacrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_shiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_trace_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_trace_curvrhs.");
  return;
}




/* DGFE Definitions */

#include <hrscc.hh>

#define config_sdg_order      5
#define config_riemann_solver hrscc::LaxFriedrichsRS<DGFE_ML_BSSN_DGFE_RHSStaticBoundary, false>

/* Export definitions */
#define ML_BSSN_DGFE_RHSStaticBoundary_sdg_grid   hrscc::GNIGrid<hrscc::GLLElement<config_sdg_order> >
#define ML_BSSN_DGFE_RHSStaticBoundary_sdg_method hrscc::SDGMethod<DGFE_ML_BSSN_DGFE_RHSStaticBoundary, ML_BSSN_DGFE_RHSStaticBoundary_sdg_grid, config_riemann_solver>

/*** Numerical scheme ***/

/* Configuration */
#define config_method ML_BSSN_DGFE_RHSStaticBoundary_sdg_method

/* Export definitions */
#define ML_BSSN_DGFE_RHSStaticBoundary_method config_method
#define ML_BSSN_DGFE_RHSStaticBoundary_solver hrscc::CLawSolver<DGFE_ML_BSSN_DGFE_RHSStaticBoundary, config_method>



class DGFE_ML_BSSN_DGFE_RHSStaticBoundary;

namespace hrscc {
  template<>
  struct traits<DGFE_ML_BSSN_DGFE_RHSStaticBoundary> {
    // All state vector variables
    enum state_t {ialpha, iA, iAt11, iAt12, iAt13, iAt22, iAt23, iAt33, iB1, iB2, iB3, ibeta1, ibeta2, ibeta3, igt11, igt12, igt13, igt22, igt23, igt33, iphi, itrK, iXt1, iXt2, iXt3, nvars};
    enum {nequations = nvars};
    enum {nexternal = 3*nvars};
    enum {nbitmasks = 0};
    static const bool pure = false;
  };
} // namespace



class DGFE_ML_BSSN_DGFE_RHSStaticBoundary: public hrscc::CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary> {
public:
  typedef hrscc::CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary> claw;
  typedef hrscc::traits<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::state_t state_t;
  typedef hrscc::traits<DGFE_ML_BSSN_DGFE_RHSStaticBoundary> variables_t;
  static const int nvars = variables_t::nvars;
  
  DGFE_ML_BSSN_DGFE_RHSStaticBoundary();
  
  inline void prim_to_all(hrscc::Observer<claw> & observer) const
  {
  }
  
  template<hrscc::policy::direction_t dir>
  inline void fluxes(hrscc::Observer<claw> & observer) const
  {
    
    CCTK_REAL fluxalphaL;
    CCTK_REAL fluxAL;
    CCTK_REAL fluxAt11L;
    CCTK_REAL fluxAt12L;
    CCTK_REAL fluxAt13L;
    CCTK_REAL fluxAt22L;
    CCTK_REAL fluxAt23L;
    CCTK_REAL fluxAt33L;
    CCTK_REAL fluxB1L;
    CCTK_REAL fluxB2L;
    CCTK_REAL fluxB3L;
    CCTK_REAL fluxbeta1L;
    CCTK_REAL fluxbeta2L;
    CCTK_REAL fluxbeta3L;
    CCTK_REAL fluxgt11L;
    CCTK_REAL fluxgt12L;
    CCTK_REAL fluxgt13L;
    CCTK_REAL fluxgt22L;
    CCTK_REAL fluxgt23L;
    CCTK_REAL fluxgt33L;
    CCTK_REAL fluxphiL;
    CCTK_REAL fluxtrKL;
    CCTK_REAL fluxXt1L;
    CCTK_REAL fluxXt2L;
    CCTK_REAL fluxXt3L;
    
    switch (dir) {
    case hrscc::policy::x: {
      fluxalphaL = observer.field[variables_t::ialpha + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAL = observer.field[variables_t::iA + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt11L = observer.field[variables_t::iAt11 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt12L = observer.field[variables_t::iAt12 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt13L = observer.field[variables_t::iAt13 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt22L = observer.field[variables_t::iAt22 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt23L = observer.field[variables_t::iAt23 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt33L = observer.field[variables_t::iAt33 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxB1L = observer.field[variables_t::iB1 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxB2L = observer.field[variables_t::iB2 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxB3L = observer.field[variables_t::iB3 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxbeta1L = observer.field[variables_t::ibeta1 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxbeta2L = observer.field[variables_t::ibeta2 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxbeta3L = observer.field[variables_t::ibeta3 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt11L = observer.field[variables_t::igt11 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt12L = observer.field[variables_t::igt12 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt13L = observer.field[variables_t::igt13 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt22L = observer.field[variables_t::igt22 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt23L = observer.field[variables_t::igt23 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt33L = observer.field[variables_t::igt33 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxphiL = observer.field[variables_t::iphi + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxtrKL = observer.field[variables_t::itrK + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxXt1L = observer.field[variables_t::iXt1 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxXt2L = observer.field[variables_t::iXt2 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxXt3L = observer.field[variables_t::iXt3 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      break;
    }
    case hrscc::policy::y: {
      fluxalphaL = observer.field[variables_t::ialpha + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAL = observer.field[variables_t::iA + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt11L = observer.field[variables_t::iAt11 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt12L = observer.field[variables_t::iAt12 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt13L = observer.field[variables_t::iAt13 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt22L = observer.field[variables_t::iAt22 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt23L = observer.field[variables_t::iAt23 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt33L = observer.field[variables_t::iAt33 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxB1L = observer.field[variables_t::iB1 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxB2L = observer.field[variables_t::iB2 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxB3L = observer.field[variables_t::iB3 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxbeta1L = observer.field[variables_t::ibeta1 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxbeta2L = observer.field[variables_t::ibeta2 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxbeta3L = observer.field[variables_t::ibeta3 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt11L = observer.field[variables_t::igt11 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt12L = observer.field[variables_t::igt12 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt13L = observer.field[variables_t::igt13 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt22L = observer.field[variables_t::igt22 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt23L = observer.field[variables_t::igt23 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt33L = observer.field[variables_t::igt33 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxphiL = observer.field[variables_t::iphi + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxtrKL = observer.field[variables_t::itrK + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxXt1L = observer.field[variables_t::iXt1 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxXt2L = observer.field[variables_t::iXt2 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxXt3L = observer.field[variables_t::iXt3 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      break;
    }
    case hrscc::policy::z: {
      fluxalphaL = observer.field[variables_t::ialpha + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAL = observer.field[variables_t::iA + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt11L = observer.field[variables_t::iAt11 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt12L = observer.field[variables_t::iAt12 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt13L = observer.field[variables_t::iAt13 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt22L = observer.field[variables_t::iAt22 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt23L = observer.field[variables_t::iAt23 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxAt33L = observer.field[variables_t::iAt33 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxB1L = observer.field[variables_t::iB1 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxB2L = observer.field[variables_t::iB2 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxB3L = observer.field[variables_t::iB3 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxbeta1L = observer.field[variables_t::ibeta1 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxbeta2L = observer.field[variables_t::ibeta2 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxbeta3L = observer.field[variables_t::ibeta3 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt11L = observer.field[variables_t::igt11 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt12L = observer.field[variables_t::igt12 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt13L = observer.field[variables_t::igt13 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt22L = observer.field[variables_t::igt22 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt23L = observer.field[variables_t::igt23 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxgt33L = observer.field[variables_t::igt33 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxphiL = observer.field[variables_t::iphi + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxtrKL = observer.field[variables_t::itrK + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxXt1L = observer.field[variables_t::iXt1 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxXt2L = observer.field[variables_t::iXt2 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      fluxXt3L = observer.field[variables_t::iXt3 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars];
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
    }
    
    observer.flux[dir][variables_t::ialpha] = fluxalphaL;
    observer.flux[dir][variables_t::iA] = fluxAL;
    observer.flux[dir][variables_t::iAt11] = fluxAt11L;
    observer.flux[dir][variables_t::iAt12] = fluxAt12L;
    observer.flux[dir][variables_t::iAt13] = fluxAt13L;
    observer.flux[dir][variables_t::iAt22] = fluxAt22L;
    observer.flux[dir][variables_t::iAt23] = fluxAt23L;
    observer.flux[dir][variables_t::iAt33] = fluxAt33L;
    observer.flux[dir][variables_t::iB1] = fluxB1L;
    observer.flux[dir][variables_t::iB2] = fluxB2L;
    observer.flux[dir][variables_t::iB3] = fluxB3L;
    observer.flux[dir][variables_t::ibeta1] = fluxbeta1L;
    observer.flux[dir][variables_t::ibeta2] = fluxbeta2L;
    observer.flux[dir][variables_t::ibeta3] = fluxbeta3L;
    observer.flux[dir][variables_t::igt11] = fluxgt11L;
    observer.flux[dir][variables_t::igt12] = fluxgt12L;
    observer.flux[dir][variables_t::igt13] = fluxgt13L;
    observer.flux[dir][variables_t::igt22] = fluxgt22L;
    observer.flux[dir][variables_t::igt23] = fluxgt23L;
    observer.flux[dir][variables_t::igt33] = fluxgt33L;
    observer.flux[dir][variables_t::iphi] = fluxphiL;
    observer.flux[dir][variables_t::itrK] = fluxtrKL;
    observer.flux[dir][variables_t::iXt1] = fluxXt1L;
    observer.flux[dir][variables_t::iXt2] = fluxXt2L;
    observer.flux[dir][variables_t::iXt3] = fluxXt3L;
  }
  
  template<hrscc::policy::direction_t dir>
  inline void eigenvalues(hrscc::Observer<claw> & observer) const
  {
    assert(0);
  }
  
  template<hrscc::policy::direction_t dir>
  inline void eig(hrscc::Observer<claw> & observer) const
  {
    assert(0);
  }
};



namespace hrscc {
  template<> int CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = {};
  template<> int CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = {};
  template<> int CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = {};
  template<> int CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[3*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = {};
  template<> int CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::bitmask_idx[0] = {};
} // namespace hrscc



namespace {
  int varindex(const char* const varname)
  {
    const int vi = CCTK_VarIndex(varname);
    if (vi<0) CCTK_ERROR("Internal error");
    return vi;
  }
}

DGFE_ML_BSSN_DGFE_RHSStaticBoundary::DGFE_ML_BSSN_DGFE_RHSStaticBoundary()
{
  using namespace hrscc;

  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::ialpha] = varindex(CCTK_THORNSTRING "::alpha");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::iA] = varindex(CCTK_THORNSTRING "::A");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::iAt11] = varindex(CCTK_THORNSTRING "::At11");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::iAt12] = varindex(CCTK_THORNSTRING "::At12");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::iAt13] = varindex(CCTK_THORNSTRING "::At13");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::iAt22] = varindex(CCTK_THORNSTRING "::At22");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::iAt23] = varindex(CCTK_THORNSTRING "::At23");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::iAt33] = varindex(CCTK_THORNSTRING "::At33");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::iB1] = varindex(CCTK_THORNSTRING "::B1");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::iB2] = varindex(CCTK_THORNSTRING "::B2");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::iB3] = varindex(CCTK_THORNSTRING "::B3");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::ibeta1] = varindex(CCTK_THORNSTRING "::beta1");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::ibeta2] = varindex(CCTK_THORNSTRING "::beta2");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::ibeta3] = varindex(CCTK_THORNSTRING "::beta3");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::igt11] = varindex(CCTK_THORNSTRING "::gt11");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::igt12] = varindex(CCTK_THORNSTRING "::gt12");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::igt13] = varindex(CCTK_THORNSTRING "::gt13");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::igt22] = varindex(CCTK_THORNSTRING "::gt22");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::igt23] = varindex(CCTK_THORNSTRING "::gt23");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::igt33] = varindex(CCTK_THORNSTRING "::gt33");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::iphi] = varindex(CCTK_THORNSTRING "::phi");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::itrK] = varindex(CCTK_THORNSTRING "::trK");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::iXt1] = varindex(CCTK_THORNSTRING "::Xt1");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::iXt2] = varindex(CCTK_THORNSTRING "::Xt2");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::conserved_idx[variables_t::iXt3] = varindex(CCTK_THORNSTRING "::Xt3");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::ialpha] = varindex(CCTK_THORNSTRING "::alpha");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::iA] = varindex(CCTK_THORNSTRING "::A");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::iAt11] = varindex(CCTK_THORNSTRING "::At11");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::iAt12] = varindex(CCTK_THORNSTRING "::At12");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::iAt13] = varindex(CCTK_THORNSTRING "::At13");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::iAt22] = varindex(CCTK_THORNSTRING "::At22");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::iAt23] = varindex(CCTK_THORNSTRING "::At23");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::iAt33] = varindex(CCTK_THORNSTRING "::At33");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::iB1] = varindex(CCTK_THORNSTRING "::B1");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::iB2] = varindex(CCTK_THORNSTRING "::B2");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::iB3] = varindex(CCTK_THORNSTRING "::B3");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::ibeta1] = varindex(CCTK_THORNSTRING "::beta1");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::ibeta2] = varindex(CCTK_THORNSTRING "::beta2");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::ibeta3] = varindex(CCTK_THORNSTRING "::beta3");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::igt11] = varindex(CCTK_THORNSTRING "::gt11");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::igt12] = varindex(CCTK_THORNSTRING "::gt12");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::igt13] = varindex(CCTK_THORNSTRING "::gt13");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::igt22] = varindex(CCTK_THORNSTRING "::gt22");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::igt23] = varindex(CCTK_THORNSTRING "::gt23");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::igt33] = varindex(CCTK_THORNSTRING "::gt33");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::iphi] = varindex(CCTK_THORNSTRING "::phi");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::itrK] = varindex(CCTK_THORNSTRING "::trK");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::iXt1] = varindex(CCTK_THORNSTRING "::Xt1");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::iXt2] = varindex(CCTK_THORNSTRING "::Xt2");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::primitive_idx[variables_t::iXt3] = varindex(CCTK_THORNSTRING "::Xt3");

  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::ialpha + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxalpha1");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iA + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxA1");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt11 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt111");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt12 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt121");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt13 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt131");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt22 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt221");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt23 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt231");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt33 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt331");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iB1 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxB11");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iB2 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxB21");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iB3 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxB31");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::ibeta1 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta11");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::ibeta2 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta21");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::ibeta3 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta31");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt11 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt111");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt12 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt121");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt13 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt131");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt22 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt221");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt23 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt231");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt33 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt331");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iphi + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxphi1");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::itrK + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxtrK1");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iXt1 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxXt11");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iXt2 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxXt21");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iXt3 + 0*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxXt31");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::ialpha + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxalpha2");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iA + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxA2");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt11 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt112");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt12 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt122");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt13 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt132");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt22 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt222");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt23 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt232");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt33 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt332");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iB1 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxB12");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iB2 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxB22");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iB3 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxB32");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::ibeta1 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta12");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::ibeta2 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta22");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::ibeta3 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta32");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt11 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt112");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt12 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt122");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt13 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt132");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt22 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt222");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt23 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt232");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt33 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt332");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iphi + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxphi2");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::itrK + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxtrK2");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iXt1 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxXt12");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iXt2 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxXt22");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iXt3 + 1*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxXt32");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::ialpha + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxalpha3");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iA + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxA3");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt11 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt113");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt12 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt123");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt13 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt133");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt22 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt223");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt23 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt233");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iAt33 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxAt333");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iB1 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxB13");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iB2 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxB23");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iB3 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxB33");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::ibeta1 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta13");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::ibeta2 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta23");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::ibeta3 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta33");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt11 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt113");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt12 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt123");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt13 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt133");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt22 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt223");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt23 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt233");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::igt33 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxgt333");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iphi + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxphi3");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::itrK + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxtrK3");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iXt1 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxXt13");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iXt2 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxXt23");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::field_idx[variables_t::iXt3 + 2*DGFE_ML_BSSN_DGFE_RHSStaticBoundary::nvars] = varindex(CCTK_THORNSTRING "::fluxXt33");

  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::ialpha] = varindex(CCTK_THORNSTRING "::alpharhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::iA] = varindex(CCTK_THORNSTRING "::Arhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::iAt11] = varindex(CCTK_THORNSTRING "::At11rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::iAt12] = varindex(CCTK_THORNSTRING "::At12rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::iAt13] = varindex(CCTK_THORNSTRING "::At13rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::iAt22] = varindex(CCTK_THORNSTRING "::At22rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::iAt23] = varindex(CCTK_THORNSTRING "::At23rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::iAt33] = varindex(CCTK_THORNSTRING "::At33rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::iB1] = varindex(CCTK_THORNSTRING "::B1rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::iB2] = varindex(CCTK_THORNSTRING "::B2rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::iB3] = varindex(CCTK_THORNSTRING "::B3rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::ibeta1] = varindex(CCTK_THORNSTRING "::beta1rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::ibeta2] = varindex(CCTK_THORNSTRING "::beta2rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::ibeta3] = varindex(CCTK_THORNSTRING "::beta3rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::igt11] = varindex(CCTK_THORNSTRING "::gt11rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::igt12] = varindex(CCTK_THORNSTRING "::gt12rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::igt13] = varindex(CCTK_THORNSTRING "::gt13rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::igt22] = varindex(CCTK_THORNSTRING "::gt22rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::igt23] = varindex(CCTK_THORNSTRING "::gt23rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::igt33] = varindex(CCTK_THORNSTRING "::gt33rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::iphi] = varindex(CCTK_THORNSTRING "::phirhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::itrK] = varindex(CCTK_THORNSTRING "::trKrhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::iXt1] = varindex(CCTK_THORNSTRING "::Xt1rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::iXt2] = varindex(CCTK_THORNSTRING "::Xt2rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHSStaticBoundary>::rhs_idx[variables_t::iXt3] = varindex(CCTK_THORNSTRING "::Xt3rhs");
}



/* A solver, DGFE's equivalent of cctkGH */
static ML_BSSN_DGFE_RHSStaticBoundary_solver *solver = NULL;



/* Call the pointwise DGFE derivative operator */
#undef PDstandardNth1
#undef PDstandardNth2
#undef PDstandardNth3
#define PDstandardNth1(u) (solver->wdiff<hrscc::policy::x>(&(u)[-index], i,j,k))
#define PDstandardNth2(u) (solver->wdiff<hrscc::policy::y>(&(u)[-index], i,j,k))
#define PDstandardNth3(u) (solver->wdiff<hrscc::policy::z>(&(u)[-index], i,j,k))



static void ML_BSSN_DGFE_RHSStaticBoundary_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Include user-supplied include files */
  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const CCTK_REAL_VEC dx CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(0));
  const CCTK_REAL_VEC dy CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(1));
  const CCTK_REAL_VEC dz CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(2));
  const CCTK_REAL_VEC dt CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_TIME);
  const CCTK_REAL_VEC t CCTK_ATTRIBUTE_UNUSED = ToReal(cctk_time);
  const CCTK_REAL_VEC dxi CCTK_ATTRIBUTE_UNUSED = INV(dx);
  const CCTK_REAL_VEC dyi CCTK_ATTRIBUTE_UNUSED = INV(dy);
  const CCTK_REAL_VEC dzi CCTK_ATTRIBUTE_UNUSED = INV(dz);
  const CCTK_REAL_VEC khalf CCTK_ATTRIBUTE_UNUSED = ToReal(0.5);
  const CCTK_REAL_VEC kthird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(0.333333333333333333333333333333);
  const CCTK_REAL_VEC ktwothird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(0.666666666666666666666666666667);
  const CCTK_REAL_VEC kfourthird CCTK_ATTRIBUTE_UNUSED = 
    ToReal(1.33333333333333333333333333333);
  const CCTK_REAL_VEC hdxi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dxi,ToReal(0.5));
  const CCTK_REAL_VEC hdyi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dyi,ToReal(0.5));
  const CCTK_REAL_VEC hdzi CCTK_ATTRIBUTE_UNUSED = 
    kmul(dzi,ToReal(0.5));
  /* Initialize predefined quantities */
  const CCTK_REAL_VEC p1o12dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dx);
  const CCTK_REAL_VEC p1o12dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dy);
  const CCTK_REAL_VEC p1o12dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dz);
  const CCTK_REAL_VEC p1o144dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dy,dx));
  const CCTK_REAL_VEC p1o144dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dz,dx));
  const CCTK_REAL_VEC p1o144dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dz,dy));
  const CCTK_REAL_VEC p1o24dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dx);
  const CCTK_REAL_VEC p1o24dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dy);
  const CCTK_REAL_VEC p1o24dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dz);
  const CCTK_REAL_VEC p1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dx);
  const CCTK_REAL_VEC p1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dy);
  const CCTK_REAL_VEC p1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dz);
  const CCTK_REAL_VEC p1o4dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dx);
  const CCTK_REAL_VEC p1o4dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dy,dx));
  const CCTK_REAL_VEC p1o4dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dz,dx));
  const CCTK_REAL_VEC p1o4dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dy);
  const CCTK_REAL_VEC p1o4dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dz,dy));
  const CCTK_REAL_VEC p1o4dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dz);
  const CCTK_REAL_VEC p1o64dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dx);
  const CCTK_REAL_VEC p1o64dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dy);
  const CCTK_REAL_VEC p1o64dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dz);
  const CCTK_REAL_VEC p1odx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dx);
  const CCTK_REAL_VEC p1odx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dx,dx));
  const CCTK_REAL_VEC p1ody CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dy);
  const CCTK_REAL_VEC p1ody2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dy,dy));
  const CCTK_REAL_VEC p1odz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dz);
  const CCTK_REAL_VEC p1odz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dz,dz));
  const CCTK_REAL_VEC pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dx,dx));
  const CCTK_REAL_VEC pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dy,dy));
  const CCTK_REAL_VEC pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dz,dz));
  const CCTK_REAL_VEC pm1o16dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dx);
  const CCTK_REAL_VEC pm1o16dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dy);
  const CCTK_REAL_VEC pm1o16dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dz);
  const CCTK_REAL_VEC pm1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dx);
  const CCTK_REAL_VEC pm1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dy);
  const CCTK_REAL_VEC pm1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dz);
  const CCTK_REAL_VEC pm1o4dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dx);
  const CCTK_REAL_VEC pm1o4dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dy);
  const CCTK_REAL_VEC pm1o4dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dz);
  /* Jacobian variable pointers */
  const bool use_jacobian1 = (!CCTK_IsFunctionAliased("MultiPatch_GetMap") || MultiPatch_GetMap(cctkGH) != jacobian_identity_map)
                        && strlen(jacobian_group) > 0;
  const bool use_jacobian = assume_use_jacobian>=0 ? assume_use_jacobian : use_jacobian1;
  const bool usejacobian CCTK_ATTRIBUTE_UNUSED = use_jacobian;
  if (use_jacobian && (strlen(jacobian_derivative_group) == 0))
  {
    CCTK_WARN(1, "GenericFD::jacobian_group and GenericFD::jacobian_derivative_group must both be set to valid group names");
  }
  
  const CCTK_REAL* restrict jacobian_ptrs[9];
  if (use_jacobian) GroupDataPointers(cctkGH, jacobian_group,
                                                9, jacobian_ptrs);
  
  const CCTK_REAL* restrict const J11 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[0] : 0;
  const CCTK_REAL* restrict const J12 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[1] : 0;
  const CCTK_REAL* restrict const J13 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[2] : 0;
  const CCTK_REAL* restrict const J21 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[3] : 0;
  const CCTK_REAL* restrict const J22 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[4] : 0;
  const CCTK_REAL* restrict const J23 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[5] : 0;
  const CCTK_REAL* restrict const J31 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[6] : 0;
  const CCTK_REAL* restrict const J32 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[7] : 0;
  const CCTK_REAL* restrict const J33 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[8] : 0;
  
  const CCTK_REAL* restrict jacobian_determinant_ptrs[1] CCTK_ATTRIBUTE_UNUSED;
  if (use_jacobian && strlen(jacobian_determinant_group) > 0) GroupDataPointers(cctkGH, jacobian_determinant_group,
                                                1, jacobian_determinant_ptrs);
  
  const CCTK_REAL* restrict const detJ CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[0] : 0;
  
  const CCTK_REAL* restrict jacobian_inverse_ptrs[9] CCTK_ATTRIBUTE_UNUSED;
  if (use_jacobian && strlen(jacobian_inverse_group) > 0) GroupDataPointers(cctkGH, jacobian_inverse_group,
                                                9, jacobian_inverse_ptrs);
  
  const CCTK_REAL* restrict const iJ11 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[0] : 0;
  const CCTK_REAL* restrict const iJ12 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[1] : 0;
  const CCTK_REAL* restrict const iJ13 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[2] : 0;
  const CCTK_REAL* restrict const iJ21 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[3] : 0;
  const CCTK_REAL* restrict const iJ22 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[4] : 0;
  const CCTK_REAL* restrict const iJ23 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[5] : 0;
  const CCTK_REAL* restrict const iJ31 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[6] : 0;
  const CCTK_REAL* restrict const iJ32 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[7] : 0;
  const CCTK_REAL* restrict const iJ33 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[8] : 0;
  
  const CCTK_REAL* restrict jacobian_derivative_ptrs[18] CCTK_ATTRIBUTE_UNUSED;
  if (use_jacobian) GroupDataPointers(cctkGH, jacobian_derivative_group,
                                      18, jacobian_derivative_ptrs);
  
  const CCTK_REAL* restrict const dJ111 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[0] : 0;
  const CCTK_REAL* restrict const dJ112 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[1] : 0;
  const CCTK_REAL* restrict const dJ113 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[2] : 0;
  const CCTK_REAL* restrict const dJ122 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[3] : 0;
  const CCTK_REAL* restrict const dJ123 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[4] : 0;
  const CCTK_REAL* restrict const dJ133 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[5] : 0;
  const CCTK_REAL* restrict const dJ211 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[6] : 0;
  const CCTK_REAL* restrict const dJ212 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[7] : 0;
  const CCTK_REAL* restrict const dJ213 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[8] : 0;
  const CCTK_REAL* restrict const dJ222 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[9] : 0;
  const CCTK_REAL* restrict const dJ223 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[10] : 0;
  const CCTK_REAL* restrict const dJ233 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[11] : 0;
  const CCTK_REAL* restrict const dJ311 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[12] : 0;
  const CCTK_REAL* restrict const dJ312 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[13] : 0;
  const CCTK_REAL* restrict const dJ313 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[14] : 0;
  const CCTK_REAL* restrict const dJ322 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[15] : 0;
  const CCTK_REAL* restrict const dJ323 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[16] : 0;
  const CCTK_REAL* restrict const dJ333 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[17] : 0;
  /* Assign local copies of arrays functions */
  
  
  /* Calculate temporaries and arrays functions */
  /* Copy local copies back to grid functions */
  /* Loop over the grid points */
  const int imin0=imin[0];
  const int imin1=imin[1];
  const int imin2=imin[2];
  const int imax0=imax[0];
  const int imax1=imax[1];
  const int imax2=imax[2];
  #pragma omp parallel
  CCTK_LOOP3STR(ML_BSSN_DGFE_RHSStaticBoundary,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    
    switch (fdOrder)
    {
      case 2:
      {
        break;
      }
      
      case 4:
      {
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL_VEC phirhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC gt11rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC gt12rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC gt13rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC gt22rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC gt23rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC gt33rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC trKrhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At11rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At12rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At13rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At22rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At23rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At33rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Xt1rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Xt2rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Xt3rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC alpharhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC ArhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC beta1rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC beta2rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC beta3rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC B1rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC B2rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC B3rhsL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(alpharhs[index],alpharhsL);
    vec_store_nta_partial(Arhs[index],ArhsL);
    vec_store_nta_partial(At11rhs[index],At11rhsL);
    vec_store_nta_partial(At12rhs[index],At12rhsL);
    vec_store_nta_partial(At13rhs[index],At13rhsL);
    vec_store_nta_partial(At22rhs[index],At22rhsL);
    vec_store_nta_partial(At23rhs[index],At23rhsL);
    vec_store_nta_partial(At33rhs[index],At33rhsL);
    vec_store_nta_partial(B1rhs[index],B1rhsL);
    vec_store_nta_partial(B2rhs[index],B2rhsL);
    vec_store_nta_partial(B3rhs[index],B3rhsL);
    vec_store_nta_partial(beta1rhs[index],beta1rhsL);
    vec_store_nta_partial(beta2rhs[index],beta2rhsL);
    vec_store_nta_partial(beta3rhs[index],beta3rhsL);
    vec_store_nta_partial(gt11rhs[index],gt11rhsL);
    vec_store_nta_partial(gt12rhs[index],gt12rhsL);
    vec_store_nta_partial(gt13rhs[index],gt13rhsL);
    vec_store_nta_partial(gt22rhs[index],gt22rhsL);
    vec_store_nta_partial(gt23rhs[index],gt23rhsL);
    vec_store_nta_partial(gt33rhs[index],gt33rhsL);
    vec_store_nta_partial(phirhs[index],phirhsL);
    vec_store_nta_partial(trKrhs[index],trKrhsL);
    vec_store_nta_partial(Xt1rhs[index],Xt1rhsL);
    vec_store_nta_partial(Xt2rhs[index],Xt2rhsL);
    vec_store_nta_partial(Xt3rhs[index],Xt3rhsL);
  }
  CCTK_ENDLOOP3STR(ML_BSSN_DGFE_RHSStaticBoundary);
}
extern "C" void ML_BSSN_DGFE_RHSStaticBoundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_DGFE_RHSStaticBoundary_Body");
  }
  if (cctk_iteration % ML_BSSN_DGFE_RHSStaticBoundary_calc_every != ML_BSSN_DGFE_RHSStaticBoundary_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ML_BSSN_DGFE::ML_curvrhs",
    "ML_BSSN_DGFE::ML_dtlapserhs",
    "ML_BSSN_DGFE::ML_dtshiftrhs",
    "ML_BSSN_DGFE::ML_Gammarhs",
    "ML_BSSN_DGFE::ML_lapserhs",
    "ML_BSSN_DGFE::ML_log_confacrhs",
    "ML_BSSN_DGFE::ML_metricrhs",
    "ML_BSSN_DGFE::ML_shiftrhs",
    "ML_BSSN_DGFE::ML_trace_curvrhs"};
  AssertGroupStorage(cctkGH, "ML_BSSN_DGFE_RHSStaticBoundary", 9, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      break;
    }
    
    case 4:
    {
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  
  if (not solver) solver = new ML_BSSN_DGFE_RHSStaticBoundary_method(cctkGH);
  LoopOverBoundary(cctkGH, ML_BSSN_DGFE_RHSStaticBoundary_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_DGFE_RHSStaticBoundary_Body");
  }
}

} // namespace ML_BSSN_DGFE
