/*  File produced by Kranc */

#define KRANC_C

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "GenericFD.h"
#include "Differencing.h"
#include "loopcontrol.h"
#include "Kranc.hh"
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

extern "C" void ML_BSSN_DGFE_RHS1_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_DGFE_RHS1_calc_every != ML_BSSN_DGFE_RHS1_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_dtlapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_dtlapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_dtshiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_Gammarhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_log_confacrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_shiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_trace_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_trace_curvrhs.");
  return;
}



/* DGFE Definitions */

#include <hrscc.hh>

#define config_sdg_order      5
#define config_riemann_solver hrscc::LaxFriedrichsRS<DGFE_ML_BSSN_DGFE_RHS1, false>

/* Export definitions */
#define ML_BSSN_DGFE_RHS1_sdg_grid   hrscc::GNIGrid<hrscc::GLLElement<config_sdg_order> >
#define ML_BSSN_DGFE_RHS1_sdg_method hrscc::SDGMethod<DGFE_ML_BSSN_DGFE_RHS1, ML_BSSN_DGFE_RHS1_sdg_grid, config_riemann_solver>

/*** Numerical scheme ***/

/* Configuration */
#define config_method ML_BSSN_DGFE_RHS1_sdg_method

/* Export definitions */
#define ML_BSSN_DGFE_RHS1_method config_method
#define ML_BSSN_DGFE_RHS1_solver hrscc::CLawSolver<DGFE_ML_BSSN_DGFE_RHS1, config_method>



class DGFE_ML_BSSN_DGFE_RHS1;

namespace hrscc {
  template<>
  struct traits<DGFE_ML_BSSN_DGFE_RHS1> {
    // All state vector variables
    enum state_t {ialpha, iA, iB1, iB2, iB3, ibeta1, ibeta2, ibeta3, igt11, igt12, igt13, igt22, igt23, igt33, iphi, itrK, iXt1, iXt2, iXt3, nvars};
    enum {nequations = nvars};
    enum {nexternal = 3*nvars};
    enum {nbitmasks = 0};
    static const bool pure = false;
  };
} // namespace



class DGFE_ML_BSSN_DGFE_RHS1: public hrscc::CLaw<DGFE_ML_BSSN_DGFE_RHS1> {
public:
  typedef hrscc::CLaw<DGFE_ML_BSSN_DGFE_RHS1> claw;
  typedef hrscc::traits<DGFE_ML_BSSN_DGFE_RHS1>::state_t state_t;
  typedef hrscc::traits<DGFE_ML_BSSN_DGFE_RHS1> variables_t;
  static const int nvars = variables_t::nvars;
  
  DGFE_ML_BSSN_DGFE_RHS1();
  
  inline void prim_to_all(hrscc::Observer<claw> & observer) const
  {
  }
  
  template<hrscc::policy::direction_t dir>
  inline void fluxes(hrscc::Observer<claw> & observer) const
  {
    
    CCTK_REAL fluxalphaL;
    CCTK_REAL fluxAL;
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
      fluxalphaL = observer.field[variables_t::ialpha + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxAL = observer.field[variables_t::iA + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxB1L = observer.field[variables_t::iB1 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxB2L = observer.field[variables_t::iB2 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxB3L = observer.field[variables_t::iB3 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxbeta1L = observer.field[variables_t::ibeta1 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxbeta2L = observer.field[variables_t::ibeta2 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxbeta3L = observer.field[variables_t::ibeta3 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt11L = observer.field[variables_t::igt11 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt12L = observer.field[variables_t::igt12 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt13L = observer.field[variables_t::igt13 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt22L = observer.field[variables_t::igt22 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt23L = observer.field[variables_t::igt23 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt33L = observer.field[variables_t::igt33 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxphiL = observer.field[variables_t::iphi + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxtrKL = observer.field[variables_t::itrK + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxXt1L = observer.field[variables_t::iXt1 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxXt2L = observer.field[variables_t::iXt2 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxXt3L = observer.field[variables_t::iXt3 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      break;
    }
    case hrscc::policy::y: {
      fluxalphaL = observer.field[variables_t::ialpha + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxAL = observer.field[variables_t::iA + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxB1L = observer.field[variables_t::iB1 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxB2L = observer.field[variables_t::iB2 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxB3L = observer.field[variables_t::iB3 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxbeta1L = observer.field[variables_t::ibeta1 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxbeta2L = observer.field[variables_t::ibeta2 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxbeta3L = observer.field[variables_t::ibeta3 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt11L = observer.field[variables_t::igt11 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt12L = observer.field[variables_t::igt12 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt13L = observer.field[variables_t::igt13 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt22L = observer.field[variables_t::igt22 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt23L = observer.field[variables_t::igt23 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt33L = observer.field[variables_t::igt33 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxphiL = observer.field[variables_t::iphi + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxtrKL = observer.field[variables_t::itrK + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxXt1L = observer.field[variables_t::iXt1 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxXt2L = observer.field[variables_t::iXt2 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxXt3L = observer.field[variables_t::iXt3 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      break;
    }
    case hrscc::policy::z: {
      fluxalphaL = observer.field[variables_t::ialpha + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxAL = observer.field[variables_t::iA + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxB1L = observer.field[variables_t::iB1 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxB2L = observer.field[variables_t::iB2 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxB3L = observer.field[variables_t::iB3 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxbeta1L = observer.field[variables_t::ibeta1 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxbeta2L = observer.field[variables_t::ibeta2 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxbeta3L = observer.field[variables_t::ibeta3 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt11L = observer.field[variables_t::igt11 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt12L = observer.field[variables_t::igt12 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt13L = observer.field[variables_t::igt13 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt22L = observer.field[variables_t::igt22 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt23L = observer.field[variables_t::igt23 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxgt33L = observer.field[variables_t::igt33 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxphiL = observer.field[variables_t::iphi + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxtrKL = observer.field[variables_t::itrK + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxXt1L = observer.field[variables_t::iXt1 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxXt2L = observer.field[variables_t::iXt2 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      fluxXt3L = observer.field[variables_t::iXt3 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars];
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
    }
    
    observer.flux[dir][variables_t::ialpha] = fluxalphaL;
    observer.flux[dir][variables_t::iA] = fluxAL;
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
  template<> int CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[DGFE_ML_BSSN_DGFE_RHS1::nvars] = {};
  template<> int CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[DGFE_ML_BSSN_DGFE_RHS1::nvars] = {};
  template<> int CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[DGFE_ML_BSSN_DGFE_RHS1::nvars] = {};
  template<> int CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[3*DGFE_ML_BSSN_DGFE_RHS1::nvars] = {};
  template<> int CLaw<DGFE_ML_BSSN_DGFE_RHS1>::bitmask_idx[0] = {};
} // namespace hrscc



namespace {
  int varindex(const char* const varname)
  {
    const int vi = CCTK_VarIndex(varname);
    if (vi<0) CCTK_ERROR("Internal error");
    return vi;
  }
}

DGFE_ML_BSSN_DGFE_RHS1::DGFE_ML_BSSN_DGFE_RHS1()
{
  using namespace hrscc;

  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::ialpha] = varindex(CCTK_THORNSTRING "::alpha");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::iA] = varindex(CCTK_THORNSTRING "::A");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::iB1] = varindex(CCTK_THORNSTRING "::B1");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::iB2] = varindex(CCTK_THORNSTRING "::B2");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::iB3] = varindex(CCTK_THORNSTRING "::B3");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::ibeta1] = varindex(CCTK_THORNSTRING "::beta1");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::ibeta2] = varindex(CCTK_THORNSTRING "::beta2");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::ibeta3] = varindex(CCTK_THORNSTRING "::beta3");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::igt11] = varindex(CCTK_THORNSTRING "::gt11");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::igt12] = varindex(CCTK_THORNSTRING "::gt12");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::igt13] = varindex(CCTK_THORNSTRING "::gt13");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::igt22] = varindex(CCTK_THORNSTRING "::gt22");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::igt23] = varindex(CCTK_THORNSTRING "::gt23");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::igt33] = varindex(CCTK_THORNSTRING "::gt33");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::iphi] = varindex(CCTK_THORNSTRING "::phi");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::itrK] = varindex(CCTK_THORNSTRING "::trK");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::iXt1] = varindex(CCTK_THORNSTRING "::Xt1");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::iXt2] = varindex(CCTK_THORNSTRING "::Xt2");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::conserved_idx[variables_t::iXt3] = varindex(CCTK_THORNSTRING "::Xt3");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::ialpha] = varindex(CCTK_THORNSTRING "::alpha");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::iA] = varindex(CCTK_THORNSTRING "::A");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::iB1] = varindex(CCTK_THORNSTRING "::B1");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::iB2] = varindex(CCTK_THORNSTRING "::B2");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::iB3] = varindex(CCTK_THORNSTRING "::B3");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::ibeta1] = varindex(CCTK_THORNSTRING "::beta1");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::ibeta2] = varindex(CCTK_THORNSTRING "::beta2");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::ibeta3] = varindex(CCTK_THORNSTRING "::beta3");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::igt11] = varindex(CCTK_THORNSTRING "::gt11");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::igt12] = varindex(CCTK_THORNSTRING "::gt12");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::igt13] = varindex(CCTK_THORNSTRING "::gt13");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::igt22] = varindex(CCTK_THORNSTRING "::gt22");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::igt23] = varindex(CCTK_THORNSTRING "::gt23");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::igt33] = varindex(CCTK_THORNSTRING "::gt33");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::iphi] = varindex(CCTK_THORNSTRING "::phi");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::itrK] = varindex(CCTK_THORNSTRING "::trK");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::iXt1] = varindex(CCTK_THORNSTRING "::Xt1");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::iXt2] = varindex(CCTK_THORNSTRING "::Xt2");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::primitive_idx[variables_t::iXt3] = varindex(CCTK_THORNSTRING "::Xt3");

  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::ialpha + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxalpha1");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iA + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxA1");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iB1 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxB11");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iB2 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxB21");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iB3 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxB31");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::ibeta1 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta11");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::ibeta2 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta21");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::ibeta3 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta31");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt11 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt111");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt12 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt121");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt13 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt131");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt22 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt221");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt23 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt231");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt33 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt331");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iphi + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxphi1");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::itrK + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxtrK1");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iXt1 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxXt11");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iXt2 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxXt21");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iXt3 + 0*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxXt31");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::ialpha + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxalpha2");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iA + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxA2");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iB1 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxB12");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iB2 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxB22");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iB3 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxB32");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::ibeta1 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta12");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::ibeta2 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta22");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::ibeta3 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta32");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt11 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt112");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt12 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt122");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt13 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt132");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt22 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt222");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt23 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt232");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt33 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt332");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iphi + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxphi2");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::itrK + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxtrK2");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iXt1 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxXt12");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iXt2 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxXt22");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iXt3 + 1*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxXt32");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::ialpha + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxalpha3");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iA + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxA3");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iB1 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxB13");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iB2 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxB23");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iB3 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxB33");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::ibeta1 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta13");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::ibeta2 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta23");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::ibeta3 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxbeta33");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt11 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt113");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt12 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt123");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt13 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt133");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt22 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt223");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt23 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt233");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::igt33 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxgt333");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iphi + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxphi3");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::itrK + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxtrK3");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iXt1 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxXt13");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iXt2 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxXt23");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::field_idx[variables_t::iXt3 + 2*DGFE_ML_BSSN_DGFE_RHS1::nvars] = varindex(CCTK_THORNSTRING "::fluxXt33");

  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::ialpha] = varindex(CCTK_THORNSTRING "::alpharhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::iA] = varindex(CCTK_THORNSTRING "::Arhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::iB1] = varindex(CCTK_THORNSTRING "::B1rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::iB2] = varindex(CCTK_THORNSTRING "::B2rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::iB3] = varindex(CCTK_THORNSTRING "::B3rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::ibeta1] = varindex(CCTK_THORNSTRING "::beta1rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::ibeta2] = varindex(CCTK_THORNSTRING "::beta2rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::ibeta3] = varindex(CCTK_THORNSTRING "::beta3rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::igt11] = varindex(CCTK_THORNSTRING "::gt11rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::igt12] = varindex(CCTK_THORNSTRING "::gt12rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::igt13] = varindex(CCTK_THORNSTRING "::gt13rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::igt22] = varindex(CCTK_THORNSTRING "::gt22rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::igt23] = varindex(CCTK_THORNSTRING "::gt23rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::igt33] = varindex(CCTK_THORNSTRING "::gt33rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::iphi] = varindex(CCTK_THORNSTRING "::phirhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::itrK] = varindex(CCTK_THORNSTRING "::trKrhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::iXt1] = varindex(CCTK_THORNSTRING "::Xt1rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::iXt2] = varindex(CCTK_THORNSTRING "::Xt2rhs");
  CLaw<DGFE_ML_BSSN_DGFE_RHS1>::rhs_idx[variables_t::iXt3] = varindex(CCTK_THORNSTRING "::Xt3rhs");
}



/* A solver, DGFE's equivalent of cctkGH */
static ML_BSSN_DGFE_RHS1_solver *solver = NULL;



/* Call the pointwise DGFE derivative operator */
#undef PDstandardNth1
#undef PDstandardNth2
#undef PDstandardNth3
#define PDstandardNth1(u) (solver->wdiff<hrscc::policy::x>(&(u)[-index], i,j,k))
#define PDstandardNth2(u) (solver->wdiff<hrscc::policy::y>(&(u)[-index], i,j,k))
#define PDstandardNth3(u) (solver->wdiff<hrscc::policy::z>(&(u)[-index], i,j,k))




static void ML_BSSN_DGFE_RHS1_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  if (use_jacobian) GenericFD_GroupDataPointers(cctkGH, jacobian_group,
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
  if (use_jacobian && strlen(jacobian_determinant_group) > 0) GenericFD_GroupDataPointers(cctkGH, jacobian_determinant_group,
                                                1, jacobian_determinant_ptrs);
  
  const CCTK_REAL* restrict const detJ CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[0] : 0;
  
  const CCTK_REAL* restrict jacobian_inverse_ptrs[9] CCTK_ATTRIBUTE_UNUSED;
  if (use_jacobian && strlen(jacobian_inverse_group) > 0) GenericFD_GroupDataPointers(cctkGH, jacobian_inverse_group,
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
  if (use_jacobian) GenericFD_GroupDataPointers(cctkGH, jacobian_derivative_group,
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
  #pragma omp parallel // reduction(+: vec_iter_counter, vec_op_counter, vec_mem_counter)
  CCTK_LOOP3STR(ML_BSSN_DGFE_RHS1,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    // vec_iter_counter+=CCTK_REAL_VEC_SIZE;
    
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC AL CCTK_ATTRIBUTE_UNUSED = vec_load(A[index]);
    CCTK_REAL_VEC alphaL CCTK_ATTRIBUTE_UNUSED = vec_load(alpha[index]);
    CCTK_REAL_VEC At11L CCTK_ATTRIBUTE_UNUSED = vec_load(At11[index]);
    CCTK_REAL_VEC At12L CCTK_ATTRIBUTE_UNUSED = vec_load(At12[index]);
    CCTK_REAL_VEC At13L CCTK_ATTRIBUTE_UNUSED = vec_load(At13[index]);
    CCTK_REAL_VEC At22L CCTK_ATTRIBUTE_UNUSED = vec_load(At22[index]);
    CCTK_REAL_VEC At23L CCTK_ATTRIBUTE_UNUSED = vec_load(At23[index]);
    CCTK_REAL_VEC At33L CCTK_ATTRIBUTE_UNUSED = vec_load(At33[index]);
    CCTK_REAL_VEC B1L CCTK_ATTRIBUTE_UNUSED = vec_load(B1[index]);
    CCTK_REAL_VEC B2L CCTK_ATTRIBUTE_UNUSED = vec_load(B2[index]);
    CCTK_REAL_VEC B3L CCTK_ATTRIBUTE_UNUSED = vec_load(B3[index]);
    CCTK_REAL_VEC beta1L CCTK_ATTRIBUTE_UNUSED = vec_load(beta1[index]);
    CCTK_REAL_VEC beta2L CCTK_ATTRIBUTE_UNUSED = vec_load(beta2[index]);
    CCTK_REAL_VEC beta3L CCTK_ATTRIBUTE_UNUSED = vec_load(beta3[index]);
    CCTK_REAL_VEC gt11L CCTK_ATTRIBUTE_UNUSED = vec_load(gt11[index]);
    CCTK_REAL_VEC gt12L CCTK_ATTRIBUTE_UNUSED = vec_load(gt12[index]);
    CCTK_REAL_VEC gt13L CCTK_ATTRIBUTE_UNUSED = vec_load(gt13[index]);
    CCTK_REAL_VEC gt22L CCTK_ATTRIBUTE_UNUSED = vec_load(gt22[index]);
    CCTK_REAL_VEC gt23L CCTK_ATTRIBUTE_UNUSED = vec_load(gt23[index]);
    CCTK_REAL_VEC gt33L CCTK_ATTRIBUTE_UNUSED = vec_load(gt33[index]);
    CCTK_REAL_VEC phiL CCTK_ATTRIBUTE_UNUSED = vec_load(phi[index]);
    CCTK_REAL_VEC rL CCTK_ATTRIBUTE_UNUSED = vec_load(r[index]);
    CCTK_REAL_VEC trKL CCTK_ATTRIBUTE_UNUSED = vec_load(trK[index]);
    CCTK_REAL_VEC Xt1L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt2L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt3L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt3[index]);
    
    CCTK_REAL_VEC eTttL, eTtxL, eTtyL, eTtzL, eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL CCTK_ATTRIBUTE_UNUSED ;
    
    if (assume_stress_energy_state>=0 ? assume_stress_energy_state : *stress_energy_state)
    {
      eTttL = vec_load(eTtt[index]);
      eTtxL = vec_load(eTtx[index]);
      eTtyL = vec_load(eTty[index]);
      eTtzL = vec_load(eTtz[index]);
      eTxxL = vec_load(eTxx[index]);
      eTxyL = vec_load(eTxy[index]);
      eTxzL = vec_load(eTxz[index]);
      eTyyL = vec_load(eTyy[index]);
      eTyzL = vec_load(eTyz[index]);
      eTzzL = vec_load(eTzz[index]);
    }
    else
    {
      eTttL = ToReal(0.0);
      eTtxL = ToReal(0.0);
      eTtyL = ToReal(0.0);
      eTtzL = ToReal(0.0);
      eTxxL = ToReal(0.0);
      eTxyL = ToReal(0.0);
      eTxzL = ToReal(0.0);
      eTyyL = ToReal(0.0);
      eTyzL = ToReal(0.0);
      eTzzL = ToReal(0.0);
    }
    
    CCTK_REAL_VEC dJ111L, dJ112L, dJ113L, dJ122L, dJ123L, dJ133L, dJ211L, dJ212L, dJ213L, dJ222L, dJ223L, dJ233L, dJ311L, dJ312L, dJ313L, dJ322L, dJ323L, dJ333L, J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L CCTK_ATTRIBUTE_UNUSED ;
    
    if (use_jacobian)
    {
      dJ111L = vec_load(dJ111[index]);
      dJ112L = vec_load(dJ112[index]);
      dJ113L = vec_load(dJ113[index]);
      dJ122L = vec_load(dJ122[index]);
      dJ123L = vec_load(dJ123[index]);
      dJ133L = vec_load(dJ133[index]);
      dJ211L = vec_load(dJ211[index]);
      dJ212L = vec_load(dJ212[index]);
      dJ213L = vec_load(dJ213[index]);
      dJ222L = vec_load(dJ222[index]);
      dJ223L = vec_load(dJ223[index]);
      dJ233L = vec_load(dJ233[index]);
      dJ311L = vec_load(dJ311[index]);
      dJ312L = vec_load(dJ312[index]);
      dJ313L = vec_load(dJ313[index]);
      dJ322L = vec_load(dJ322[index]);
      dJ323L = vec_load(dJ323[index]);
      dJ333L = vec_load(dJ333[index]);
      J11L = vec_load(J11[index]);
      J12L = vec_load(J12[index]);
      J13L = vec_load(J13[index]);
      J21L = vec_load(J21[index]);
      J22L = vec_load(J22[index]);
      J23L = vec_load(J23[index]);
      J31L = vec_load(J31[index]);
      J32L = vec_load(J32[index]);
      J33L = vec_load(J33[index]);
    }
    
    /* Include user supplied include files */
    
    /* Precompute derivatives */
    CCTK_REAL_VEC PDstandardNth1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3trK CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDstandardNth1alpha = PDstandardNthfdOrder21(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder22(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder23(&alpha[index]);
        PDstandardNth11alpha = PDstandardNthfdOrder211(&alpha[index]);
        PDstandardNth22alpha = PDstandardNthfdOrder222(&alpha[index]);
        PDstandardNth33alpha = PDstandardNthfdOrder233(&alpha[index]);
        PDstandardNth12alpha = PDstandardNthfdOrder212(&alpha[index]);
        PDstandardNth13alpha = PDstandardNthfdOrder213(&alpha[index]);
        PDstandardNth23alpha = PDstandardNthfdOrder223(&alpha[index]);
        PDstandardNth1beta1 = PDstandardNthfdOrder21(&beta1[index]);
        PDstandardNth2beta1 = PDstandardNthfdOrder22(&beta1[index]);
        PDstandardNth3beta1 = PDstandardNthfdOrder23(&beta1[index]);
        PDstandardNth11beta1 = PDstandardNthfdOrder211(&beta1[index]);
        PDstandardNth22beta1 = PDstandardNthfdOrder222(&beta1[index]);
        PDstandardNth33beta1 = PDstandardNthfdOrder233(&beta1[index]);
        PDstandardNth12beta1 = PDstandardNthfdOrder212(&beta1[index]);
        PDstandardNth13beta1 = PDstandardNthfdOrder213(&beta1[index]);
        PDstandardNth23beta1 = PDstandardNthfdOrder223(&beta1[index]);
        PDstandardNth1beta2 = PDstandardNthfdOrder21(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder22(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder23(&beta2[index]);
        PDstandardNth11beta2 = PDstandardNthfdOrder211(&beta2[index]);
        PDstandardNth22beta2 = PDstandardNthfdOrder222(&beta2[index]);
        PDstandardNth33beta2 = PDstandardNthfdOrder233(&beta2[index]);
        PDstandardNth12beta2 = PDstandardNthfdOrder212(&beta2[index]);
        PDstandardNth13beta2 = PDstandardNthfdOrder213(&beta2[index]);
        PDstandardNth23beta2 = PDstandardNthfdOrder223(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder21(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder22(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder23(&beta3[index]);
        PDstandardNth11beta3 = PDstandardNthfdOrder211(&beta3[index]);
        PDstandardNth22beta3 = PDstandardNthfdOrder222(&beta3[index]);
        PDstandardNth33beta3 = PDstandardNthfdOrder233(&beta3[index]);
        PDstandardNth12beta3 = PDstandardNthfdOrder212(&beta3[index]);
        PDstandardNth13beta3 = PDstandardNthfdOrder213(&beta3[index]);
        PDstandardNth23beta3 = PDstandardNthfdOrder223(&beta3[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder21(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder22(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder23(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder21(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder22(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder23(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder21(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder22(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder23(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder21(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder22(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder23(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder21(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder22(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder23(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder21(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder22(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder23(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder21(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder22(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder23(&phi[index]);
        PDstandardNth1trK = PDstandardNthfdOrder21(&trK[index]);
        PDstandardNth2trK = PDstandardNthfdOrder22(&trK[index]);
        PDstandardNth3trK = PDstandardNthfdOrder23(&trK[index]);
        break;
      }
      
      case 4:
      {
        PDstandardNth1alpha = PDstandardNthfdOrder41(&alpha[index]);
        PDstandardNth2alpha = PDstandardNthfdOrder42(&alpha[index]);
        PDstandardNth3alpha = PDstandardNthfdOrder43(&alpha[index]);
        PDstandardNth11alpha = PDstandardNthfdOrder411(&alpha[index]);
        PDstandardNth22alpha = PDstandardNthfdOrder422(&alpha[index]);
        PDstandardNth33alpha = PDstandardNthfdOrder433(&alpha[index]);
        PDstandardNth12alpha = PDstandardNthfdOrder412(&alpha[index]);
        PDstandardNth13alpha = PDstandardNthfdOrder413(&alpha[index]);
        PDstandardNth23alpha = PDstandardNthfdOrder423(&alpha[index]);
        PDstandardNth1beta1 = PDstandardNthfdOrder41(&beta1[index]);
        PDstandardNth2beta1 = PDstandardNthfdOrder42(&beta1[index]);
        PDstandardNth3beta1 = PDstandardNthfdOrder43(&beta1[index]);
        PDstandardNth11beta1 = PDstandardNthfdOrder411(&beta1[index]);
        PDstandardNth22beta1 = PDstandardNthfdOrder422(&beta1[index]);
        PDstandardNth33beta1 = PDstandardNthfdOrder433(&beta1[index]);
        PDstandardNth12beta1 = PDstandardNthfdOrder412(&beta1[index]);
        PDstandardNth13beta1 = PDstandardNthfdOrder413(&beta1[index]);
        PDstandardNth23beta1 = PDstandardNthfdOrder423(&beta1[index]);
        PDstandardNth1beta2 = PDstandardNthfdOrder41(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder42(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder43(&beta2[index]);
        PDstandardNth11beta2 = PDstandardNthfdOrder411(&beta2[index]);
        PDstandardNth22beta2 = PDstandardNthfdOrder422(&beta2[index]);
        PDstandardNth33beta2 = PDstandardNthfdOrder433(&beta2[index]);
        PDstandardNth12beta2 = PDstandardNthfdOrder412(&beta2[index]);
        PDstandardNth13beta2 = PDstandardNthfdOrder413(&beta2[index]);
        PDstandardNth23beta2 = PDstandardNthfdOrder423(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder41(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder42(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder43(&beta3[index]);
        PDstandardNth11beta3 = PDstandardNthfdOrder411(&beta3[index]);
        PDstandardNth22beta3 = PDstandardNthfdOrder422(&beta3[index]);
        PDstandardNth33beta3 = PDstandardNthfdOrder433(&beta3[index]);
        PDstandardNth12beta3 = PDstandardNthfdOrder412(&beta3[index]);
        PDstandardNth13beta3 = PDstandardNthfdOrder413(&beta3[index]);
        PDstandardNth23beta3 = PDstandardNthfdOrder423(&beta3[index]);
        PDstandardNth1gt11 = PDstandardNthfdOrder41(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder42(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder43(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder41(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder42(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder43(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder41(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder42(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder43(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder41(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder42(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder43(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder41(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder42(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder43(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder41(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder42(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder43(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder41(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder42(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder43(&phi[index]);
        PDstandardNth1trK = PDstandardNthfdOrder41(&trK[index]);
        PDstandardNth2trK = PDstandardNthfdOrder42(&trK[index]);
        PDstandardNth3trK = PDstandardNthfdOrder43(&trK[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = kisgn(beta1L);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = kisgn(beta2L);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = kisgn(beta3L);
    
    CCTK_REAL_VEC JacPDstandardNth11alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth11beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth11beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth11beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3trK CCTK_ATTRIBUTE_UNUSED;
    
    if (use_jacobian)
    {
      JacPDstandardNth1alpha = 
        kmadd(J11L,PDstandardNth1alpha,kmadd(J21L,PDstandardNth2alpha,kmul(J31L,PDstandardNth3alpha)));
      
      JacPDstandardNth1beta1 = 
        kmadd(J11L,PDstandardNth1beta1,kmadd(J21L,PDstandardNth2beta1,kmul(J31L,PDstandardNth3beta1)));
      
      JacPDstandardNth1beta2 = 
        kmadd(J11L,PDstandardNth1beta2,kmadd(J21L,PDstandardNth2beta2,kmul(J31L,PDstandardNth3beta2)));
      
      JacPDstandardNth1beta3 = 
        kmadd(J11L,PDstandardNth1beta3,kmadd(J21L,PDstandardNth2beta3,kmul(J31L,PDstandardNth3beta3)));
      
      JacPDstandardNth1gt11 = 
        kmadd(J11L,PDstandardNth1gt11,kmadd(J21L,PDstandardNth2gt11,kmul(J31L,PDstandardNth3gt11)));
      
      JacPDstandardNth1gt12 = 
        kmadd(J11L,PDstandardNth1gt12,kmadd(J21L,PDstandardNth2gt12,kmul(J31L,PDstandardNth3gt12)));
      
      JacPDstandardNth1gt13 = 
        kmadd(J11L,PDstandardNth1gt13,kmadd(J21L,PDstandardNth2gt13,kmul(J31L,PDstandardNth3gt13)));
      
      JacPDstandardNth1gt22 = 
        kmadd(J11L,PDstandardNth1gt22,kmadd(J21L,PDstandardNth2gt22,kmul(J31L,PDstandardNth3gt22)));
      
      JacPDstandardNth1gt23 = 
        kmadd(J11L,PDstandardNth1gt23,kmadd(J21L,PDstandardNth2gt23,kmul(J31L,PDstandardNth3gt23)));
      
      JacPDstandardNth1gt33 = 
        kmadd(J11L,PDstandardNth1gt33,kmadd(J21L,PDstandardNth2gt33,kmul(J31L,PDstandardNth3gt33)));
      
      JacPDstandardNth1phi = 
        kmadd(J11L,PDstandardNth1phi,kmadd(J21L,PDstandardNth2phi,kmul(J31L,PDstandardNth3phi)));
      
      JacPDstandardNth1trK = 
        kmadd(J11L,PDstandardNth1trK,kmadd(J21L,PDstandardNth2trK,kmul(J31L,PDstandardNth3trK)));
      
      JacPDstandardNth2alpha = 
        kmadd(J12L,PDstandardNth1alpha,kmadd(J22L,PDstandardNth2alpha,kmul(J32L,PDstandardNth3alpha)));
      
      JacPDstandardNth2beta1 = 
        kmadd(J12L,PDstandardNth1beta1,kmadd(J22L,PDstandardNth2beta1,kmul(J32L,PDstandardNth3beta1)));
      
      JacPDstandardNth2beta2 = 
        kmadd(J12L,PDstandardNth1beta2,kmadd(J22L,PDstandardNth2beta2,kmul(J32L,PDstandardNth3beta2)));
      
      JacPDstandardNth2beta3 = 
        kmadd(J12L,PDstandardNth1beta3,kmadd(J22L,PDstandardNth2beta3,kmul(J32L,PDstandardNth3beta3)));
      
      JacPDstandardNth2gt11 = 
        kmadd(J12L,PDstandardNth1gt11,kmadd(J22L,PDstandardNth2gt11,kmul(J32L,PDstandardNth3gt11)));
      
      JacPDstandardNth2gt12 = 
        kmadd(J12L,PDstandardNth1gt12,kmadd(J22L,PDstandardNth2gt12,kmul(J32L,PDstandardNth3gt12)));
      
      JacPDstandardNth2gt13 = 
        kmadd(J12L,PDstandardNth1gt13,kmadd(J22L,PDstandardNth2gt13,kmul(J32L,PDstandardNth3gt13)));
      
      JacPDstandardNth2gt22 = 
        kmadd(J12L,PDstandardNth1gt22,kmadd(J22L,PDstandardNth2gt22,kmul(J32L,PDstandardNth3gt22)));
      
      JacPDstandardNth2gt23 = 
        kmadd(J12L,PDstandardNth1gt23,kmadd(J22L,PDstandardNth2gt23,kmul(J32L,PDstandardNth3gt23)));
      
      JacPDstandardNth2gt33 = 
        kmadd(J12L,PDstandardNth1gt33,kmadd(J22L,PDstandardNth2gt33,kmul(J32L,PDstandardNth3gt33)));
      
      JacPDstandardNth2phi = 
        kmadd(J12L,PDstandardNth1phi,kmadd(J22L,PDstandardNth2phi,kmul(J32L,PDstandardNth3phi)));
      
      JacPDstandardNth2trK = 
        kmadd(J12L,PDstandardNth1trK,kmadd(J22L,PDstandardNth2trK,kmul(J32L,PDstandardNth3trK)));
      
      JacPDstandardNth3alpha = 
        kmadd(J13L,PDstandardNth1alpha,kmadd(J23L,PDstandardNth2alpha,kmul(J33L,PDstandardNth3alpha)));
      
      JacPDstandardNth3beta1 = 
        kmadd(J13L,PDstandardNth1beta1,kmadd(J23L,PDstandardNth2beta1,kmul(J33L,PDstandardNth3beta1)));
      
      JacPDstandardNth3beta2 = 
        kmadd(J13L,PDstandardNth1beta2,kmadd(J23L,PDstandardNth2beta2,kmul(J33L,PDstandardNth3beta2)));
      
      JacPDstandardNth3beta3 = 
        kmadd(J13L,PDstandardNth1beta3,kmadd(J23L,PDstandardNth2beta3,kmul(J33L,PDstandardNth3beta3)));
      
      JacPDstandardNth3gt11 = 
        kmadd(J13L,PDstandardNth1gt11,kmadd(J23L,PDstandardNth2gt11,kmul(J33L,PDstandardNth3gt11)));
      
      JacPDstandardNth3gt12 = 
        kmadd(J13L,PDstandardNth1gt12,kmadd(J23L,PDstandardNth2gt12,kmul(J33L,PDstandardNth3gt12)));
      
      JacPDstandardNth3gt13 = 
        kmadd(J13L,PDstandardNth1gt13,kmadd(J23L,PDstandardNth2gt13,kmul(J33L,PDstandardNth3gt13)));
      
      JacPDstandardNth3gt22 = 
        kmadd(J13L,PDstandardNth1gt22,kmadd(J23L,PDstandardNth2gt22,kmul(J33L,PDstandardNth3gt22)));
      
      JacPDstandardNth3gt23 = 
        kmadd(J13L,PDstandardNth1gt23,kmadd(J23L,PDstandardNth2gt23,kmul(J33L,PDstandardNth3gt23)));
      
      JacPDstandardNth3gt33 = 
        kmadd(J13L,PDstandardNth1gt33,kmadd(J23L,PDstandardNth2gt33,kmul(J33L,PDstandardNth3gt33)));
      
      JacPDstandardNth3phi = 
        kmadd(J13L,PDstandardNth1phi,kmadd(J23L,PDstandardNth2phi,kmul(J33L,PDstandardNth3phi)));
      
      JacPDstandardNth3trK = 
        kmadd(J13L,PDstandardNth1trK,kmadd(J23L,PDstandardNth2trK,kmul(J33L,PDstandardNth3trK)));
      
      JacPDstandardNth11alpha = 
        kmadd(dJ111L,PDstandardNth1alpha,kmadd(dJ211L,PDstandardNth2alpha,kmadd(dJ311L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,kmul(J11L,J11L),kmadd(PDstandardNth22alpha,kmul(J21L,J21L),kmadd(PDstandardNth33alpha,kmul(J31L,J31L),kmadd(J11L,kmul(J21L,kmul(PDstandardNth12alpha,ToReal(2))),kmadd(J11L,kmul(J31L,kmul(PDstandardNth13alpha,ToReal(2))),kmul(J21L,kmul(J31L,kmul(PDstandardNth23alpha,ToReal(2))))))))))));
      
      JacPDstandardNth11beta1 = 
        kmadd(dJ111L,PDstandardNth1beta1,kmadd(dJ211L,PDstandardNth2beta1,kmadd(dJ311L,PDstandardNth3beta1,kmadd(PDstandardNth11beta1,kmul(J11L,J11L),kmadd(PDstandardNth22beta1,kmul(J21L,J21L),kmadd(PDstandardNth33beta1,kmul(J31L,J31L),kmadd(J11L,kmul(J21L,kmul(PDstandardNth12beta1,ToReal(2))),kmadd(J11L,kmul(J31L,kmul(PDstandardNth13beta1,ToReal(2))),kmul(J21L,kmul(J31L,kmul(PDstandardNth23beta1,ToReal(2))))))))))));
      
      JacPDstandardNth11beta2 = 
        kmadd(dJ111L,PDstandardNth1beta2,kmadd(dJ211L,PDstandardNth2beta2,kmadd(dJ311L,PDstandardNth3beta2,kmadd(PDstandardNth11beta2,kmul(J11L,J11L),kmadd(PDstandardNth22beta2,kmul(J21L,J21L),kmadd(PDstandardNth33beta2,kmul(J31L,J31L),kmadd(J11L,kmul(J21L,kmul(PDstandardNth12beta2,ToReal(2))),kmadd(J11L,kmul(J31L,kmul(PDstandardNth13beta2,ToReal(2))),kmul(J21L,kmul(J31L,kmul(PDstandardNth23beta2,ToReal(2))))))))))));
      
      JacPDstandardNth11beta3 = 
        kmadd(dJ111L,PDstandardNth1beta3,kmadd(dJ211L,PDstandardNth2beta3,kmadd(dJ311L,PDstandardNth3beta3,kmadd(PDstandardNth11beta3,kmul(J11L,J11L),kmadd(PDstandardNth22beta3,kmul(J21L,J21L),kmadd(PDstandardNth33beta3,kmul(J31L,J31L),kmadd(J11L,kmul(J21L,kmul(PDstandardNth12beta3,ToReal(2))),kmadd(J11L,kmul(J31L,kmul(PDstandardNth13beta3,ToReal(2))),kmul(J21L,kmul(J31L,kmul(PDstandardNth23beta3,ToReal(2))))))))))));
      
      JacPDstandardNth22alpha = 
        kmadd(dJ122L,PDstandardNth1alpha,kmadd(dJ222L,PDstandardNth2alpha,kmadd(dJ322L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,kmul(J12L,J12L),kmadd(PDstandardNth22alpha,kmul(J22L,J22L),kmadd(PDstandardNth33alpha,kmul(J32L,J32L),kmadd(J12L,kmul(J22L,kmul(PDstandardNth12alpha,ToReal(2))),kmadd(J12L,kmul(J32L,kmul(PDstandardNth13alpha,ToReal(2))),kmul(J22L,kmul(J32L,kmul(PDstandardNth23alpha,ToReal(2))))))))))));
      
      JacPDstandardNth22beta1 = 
        kmadd(dJ122L,PDstandardNth1beta1,kmadd(dJ222L,PDstandardNth2beta1,kmadd(dJ322L,PDstandardNth3beta1,kmadd(PDstandardNth11beta1,kmul(J12L,J12L),kmadd(PDstandardNth22beta1,kmul(J22L,J22L),kmadd(PDstandardNth33beta1,kmul(J32L,J32L),kmadd(J12L,kmul(J22L,kmul(PDstandardNth12beta1,ToReal(2))),kmadd(J12L,kmul(J32L,kmul(PDstandardNth13beta1,ToReal(2))),kmul(J22L,kmul(J32L,kmul(PDstandardNth23beta1,ToReal(2))))))))))));
      
      JacPDstandardNth22beta2 = 
        kmadd(dJ122L,PDstandardNth1beta2,kmadd(dJ222L,PDstandardNth2beta2,kmadd(dJ322L,PDstandardNth3beta2,kmadd(PDstandardNth11beta2,kmul(J12L,J12L),kmadd(PDstandardNth22beta2,kmul(J22L,J22L),kmadd(PDstandardNth33beta2,kmul(J32L,J32L),kmadd(J12L,kmul(J22L,kmul(PDstandardNth12beta2,ToReal(2))),kmadd(J12L,kmul(J32L,kmul(PDstandardNth13beta2,ToReal(2))),kmul(J22L,kmul(J32L,kmul(PDstandardNth23beta2,ToReal(2))))))))))));
      
      JacPDstandardNth22beta3 = 
        kmadd(dJ122L,PDstandardNth1beta3,kmadd(dJ222L,PDstandardNth2beta3,kmadd(dJ322L,PDstandardNth3beta3,kmadd(PDstandardNth11beta3,kmul(J12L,J12L),kmadd(PDstandardNth22beta3,kmul(J22L,J22L),kmadd(PDstandardNth33beta3,kmul(J32L,J32L),kmadd(J12L,kmul(J22L,kmul(PDstandardNth12beta3,ToReal(2))),kmadd(J12L,kmul(J32L,kmul(PDstandardNth13beta3,ToReal(2))),kmul(J22L,kmul(J32L,kmul(PDstandardNth23beta3,ToReal(2))))))))))));
      
      JacPDstandardNth33alpha = 
        kmadd(dJ133L,PDstandardNth1alpha,kmadd(dJ233L,PDstandardNth2alpha,kmadd(dJ333L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,kmul(J13L,J13L),kmadd(PDstandardNth22alpha,kmul(J23L,J23L),kmadd(PDstandardNth33alpha,kmul(J33L,J33L),kmadd(J13L,kmul(J23L,kmul(PDstandardNth12alpha,ToReal(2))),kmadd(J13L,kmul(J33L,kmul(PDstandardNth13alpha,ToReal(2))),kmul(J23L,kmul(J33L,kmul(PDstandardNth23alpha,ToReal(2))))))))))));
      
      JacPDstandardNth33beta1 = 
        kmadd(dJ133L,PDstandardNth1beta1,kmadd(dJ233L,PDstandardNth2beta1,kmadd(dJ333L,PDstandardNth3beta1,kmadd(PDstandardNth11beta1,kmul(J13L,J13L),kmadd(PDstandardNth22beta1,kmul(J23L,J23L),kmadd(PDstandardNth33beta1,kmul(J33L,J33L),kmadd(J13L,kmul(J23L,kmul(PDstandardNth12beta1,ToReal(2))),kmadd(J13L,kmul(J33L,kmul(PDstandardNth13beta1,ToReal(2))),kmul(J23L,kmul(J33L,kmul(PDstandardNth23beta1,ToReal(2))))))))))));
      
      JacPDstandardNth33beta2 = 
        kmadd(dJ133L,PDstandardNth1beta2,kmadd(dJ233L,PDstandardNth2beta2,kmadd(dJ333L,PDstandardNth3beta2,kmadd(PDstandardNth11beta2,kmul(J13L,J13L),kmadd(PDstandardNth22beta2,kmul(J23L,J23L),kmadd(PDstandardNth33beta2,kmul(J33L,J33L),kmadd(J13L,kmul(J23L,kmul(PDstandardNth12beta2,ToReal(2))),kmadd(J13L,kmul(J33L,kmul(PDstandardNth13beta2,ToReal(2))),kmul(J23L,kmul(J33L,kmul(PDstandardNth23beta2,ToReal(2))))))))))));
      
      JacPDstandardNth33beta3 = 
        kmadd(dJ133L,PDstandardNth1beta3,kmadd(dJ233L,PDstandardNth2beta3,kmadd(dJ333L,PDstandardNth3beta3,kmadd(PDstandardNth11beta3,kmul(J13L,J13L),kmadd(PDstandardNth22beta3,kmul(J23L,J23L),kmadd(PDstandardNth33beta3,kmul(J33L,J33L),kmadd(J13L,kmul(J23L,kmul(PDstandardNth12beta3,ToReal(2))),kmadd(J13L,kmul(J33L,kmul(PDstandardNth13beta3,ToReal(2))),kmul(J23L,kmul(J33L,kmul(PDstandardNth23beta3,ToReal(2))))))))))));
      
      JacPDstandardNth12alpha = 
        kmadd(J11L,kmul(J12L,PDstandardNth11alpha),kmadd(J12L,kmul(J21L,PDstandardNth12alpha),kmadd(J11L,kmul(J22L,PDstandardNth12alpha),kmadd(J12L,kmul(J31L,PDstandardNth13alpha),kmadd(J11L,kmul(J32L,PDstandardNth13alpha),kmadd(dJ112L,PDstandardNth1alpha,kmadd(J21L,kmul(J22L,PDstandardNth22alpha),kmadd(J22L,kmul(J31L,PDstandardNth23alpha),kmadd(J21L,kmul(J32L,PDstandardNth23alpha),kmadd(dJ212L,PDstandardNth2alpha,kmadd(J31L,kmul(J32L,PDstandardNth33alpha),kmul(dJ312L,PDstandardNth3alpha))))))))))));
      
      JacPDstandardNth12beta1 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11beta1),kmadd(J12L,kmul(J21L,PDstandardNth12beta1),kmadd(J11L,kmul(J22L,PDstandardNth12beta1),kmadd(J12L,kmul(J31L,PDstandardNth13beta1),kmadd(J11L,kmul(J32L,PDstandardNth13beta1),kmadd(dJ112L,PDstandardNth1beta1,kmadd(J21L,kmul(J22L,PDstandardNth22beta1),kmadd(J22L,kmul(J31L,PDstandardNth23beta1),kmadd(J21L,kmul(J32L,PDstandardNth23beta1),kmadd(dJ212L,PDstandardNth2beta1,kmadd(J31L,kmul(J32L,PDstandardNth33beta1),kmul(dJ312L,PDstandardNth3beta1))))))))))));
      
      JacPDstandardNth12beta2 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11beta2),kmadd(J12L,kmul(J21L,PDstandardNth12beta2),kmadd(J11L,kmul(J22L,PDstandardNth12beta2),kmadd(J12L,kmul(J31L,PDstandardNth13beta2),kmadd(J11L,kmul(J32L,PDstandardNth13beta2),kmadd(dJ112L,PDstandardNth1beta2,kmadd(J21L,kmul(J22L,PDstandardNth22beta2),kmadd(J22L,kmul(J31L,PDstandardNth23beta2),kmadd(J21L,kmul(J32L,PDstandardNth23beta2),kmadd(dJ212L,PDstandardNth2beta2,kmadd(J31L,kmul(J32L,PDstandardNth33beta2),kmul(dJ312L,PDstandardNth3beta2))))))))))));
      
      JacPDstandardNth12beta3 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11beta3),kmadd(J12L,kmul(J21L,PDstandardNth12beta3),kmadd(J11L,kmul(J22L,PDstandardNth12beta3),kmadd(J12L,kmul(J31L,PDstandardNth13beta3),kmadd(J11L,kmul(J32L,PDstandardNth13beta3),kmadd(dJ112L,PDstandardNth1beta3,kmadd(J21L,kmul(J22L,PDstandardNth22beta3),kmadd(J22L,kmul(J31L,PDstandardNth23beta3),kmadd(J21L,kmul(J32L,PDstandardNth23beta3),kmadd(dJ212L,PDstandardNth2beta3,kmadd(J31L,kmul(J32L,PDstandardNth33beta3),kmul(dJ312L,PDstandardNth3beta3))))))))))));
      
      JacPDstandardNth13alpha = 
        kmadd(J11L,kmul(J13L,PDstandardNth11alpha),kmadd(J13L,kmul(J21L,PDstandardNth12alpha),kmadd(J11L,kmul(J23L,PDstandardNth12alpha),kmadd(J13L,kmul(J31L,PDstandardNth13alpha),kmadd(J11L,kmul(J33L,PDstandardNth13alpha),kmadd(dJ113L,PDstandardNth1alpha,kmadd(J21L,kmul(J23L,PDstandardNth22alpha),kmadd(J23L,kmul(J31L,PDstandardNth23alpha),kmadd(J21L,kmul(J33L,PDstandardNth23alpha),kmadd(dJ213L,PDstandardNth2alpha,kmadd(J31L,kmul(J33L,PDstandardNth33alpha),kmul(dJ313L,PDstandardNth3alpha))))))))))));
      
      JacPDstandardNth13beta1 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11beta1),kmadd(J13L,kmul(J21L,PDstandardNth12beta1),kmadd(J11L,kmul(J23L,PDstandardNth12beta1),kmadd(J13L,kmul(J31L,PDstandardNth13beta1),kmadd(J11L,kmul(J33L,PDstandardNth13beta1),kmadd(dJ113L,PDstandardNth1beta1,kmadd(J21L,kmul(J23L,PDstandardNth22beta1),kmadd(J23L,kmul(J31L,PDstandardNth23beta1),kmadd(J21L,kmul(J33L,PDstandardNth23beta1),kmadd(dJ213L,PDstandardNth2beta1,kmadd(J31L,kmul(J33L,PDstandardNth33beta1),kmul(dJ313L,PDstandardNth3beta1))))))))))));
      
      JacPDstandardNth13beta2 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11beta2),kmadd(J13L,kmul(J21L,PDstandardNth12beta2),kmadd(J11L,kmul(J23L,PDstandardNth12beta2),kmadd(J13L,kmul(J31L,PDstandardNth13beta2),kmadd(J11L,kmul(J33L,PDstandardNth13beta2),kmadd(dJ113L,PDstandardNth1beta2,kmadd(J21L,kmul(J23L,PDstandardNth22beta2),kmadd(J23L,kmul(J31L,PDstandardNth23beta2),kmadd(J21L,kmul(J33L,PDstandardNth23beta2),kmadd(dJ213L,PDstandardNth2beta2,kmadd(J31L,kmul(J33L,PDstandardNth33beta2),kmul(dJ313L,PDstandardNth3beta2))))))))))));
      
      JacPDstandardNth13beta3 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11beta3),kmadd(J13L,kmul(J21L,PDstandardNth12beta3),kmadd(J11L,kmul(J23L,PDstandardNth12beta3),kmadd(J13L,kmul(J31L,PDstandardNth13beta3),kmadd(J11L,kmul(J33L,PDstandardNth13beta3),kmadd(dJ113L,PDstandardNth1beta3,kmadd(J21L,kmul(J23L,PDstandardNth22beta3),kmadd(J23L,kmul(J31L,PDstandardNth23beta3),kmadd(J21L,kmul(J33L,PDstandardNth23beta3),kmadd(dJ213L,PDstandardNth2beta3,kmadd(J31L,kmul(J33L,PDstandardNth33beta3),kmul(dJ313L,PDstandardNth3beta3))))))))))));
      
      JacPDstandardNth21alpha = 
        kmadd(J11L,kmul(J12L,PDstandardNth11alpha),kmadd(J12L,kmul(J21L,PDstandardNth12alpha),kmadd(J11L,kmul(J22L,PDstandardNth12alpha),kmadd(J12L,kmul(J31L,PDstandardNth13alpha),kmadd(J11L,kmul(J32L,PDstandardNth13alpha),kmadd(dJ112L,PDstandardNth1alpha,kmadd(J21L,kmul(J22L,PDstandardNth22alpha),kmadd(J22L,kmul(J31L,PDstandardNth23alpha),kmadd(J21L,kmul(J32L,PDstandardNth23alpha),kmadd(dJ212L,PDstandardNth2alpha,kmadd(J31L,kmul(J32L,PDstandardNth33alpha),kmul(dJ312L,PDstandardNth3alpha))))))))))));
      
      JacPDstandardNth21beta1 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11beta1),kmadd(J12L,kmul(J21L,PDstandardNth12beta1),kmadd(J11L,kmul(J22L,PDstandardNth12beta1),kmadd(J12L,kmul(J31L,PDstandardNth13beta1),kmadd(J11L,kmul(J32L,PDstandardNth13beta1),kmadd(dJ112L,PDstandardNth1beta1,kmadd(J21L,kmul(J22L,PDstandardNth22beta1),kmadd(J22L,kmul(J31L,PDstandardNth23beta1),kmadd(J21L,kmul(J32L,PDstandardNth23beta1),kmadd(dJ212L,PDstandardNth2beta1,kmadd(J31L,kmul(J32L,PDstandardNth33beta1),kmul(dJ312L,PDstandardNth3beta1))))))))))));
      
      JacPDstandardNth21beta2 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11beta2),kmadd(J12L,kmul(J21L,PDstandardNth12beta2),kmadd(J11L,kmul(J22L,PDstandardNth12beta2),kmadd(J12L,kmul(J31L,PDstandardNth13beta2),kmadd(J11L,kmul(J32L,PDstandardNth13beta2),kmadd(dJ112L,PDstandardNth1beta2,kmadd(J21L,kmul(J22L,PDstandardNth22beta2),kmadd(J22L,kmul(J31L,PDstandardNth23beta2),kmadd(J21L,kmul(J32L,PDstandardNth23beta2),kmadd(dJ212L,PDstandardNth2beta2,kmadd(J31L,kmul(J32L,PDstandardNth33beta2),kmul(dJ312L,PDstandardNth3beta2))))))))))));
      
      JacPDstandardNth21beta3 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11beta3),kmadd(J12L,kmul(J21L,PDstandardNth12beta3),kmadd(J11L,kmul(J22L,PDstandardNth12beta3),kmadd(J12L,kmul(J31L,PDstandardNth13beta3),kmadd(J11L,kmul(J32L,PDstandardNth13beta3),kmadd(dJ112L,PDstandardNth1beta3,kmadd(J21L,kmul(J22L,PDstandardNth22beta3),kmadd(J22L,kmul(J31L,PDstandardNth23beta3),kmadd(J21L,kmul(J32L,PDstandardNth23beta3),kmadd(dJ212L,PDstandardNth2beta3,kmadd(J31L,kmul(J32L,PDstandardNth33beta3),kmul(dJ312L,PDstandardNth3beta3))))))))))));
      
      JacPDstandardNth23alpha = 
        kmadd(J12L,kmul(J13L,PDstandardNth11alpha),kmadd(J13L,kmul(J22L,PDstandardNth12alpha),kmadd(J12L,kmul(J23L,PDstandardNth12alpha),kmadd(J13L,kmul(J32L,PDstandardNth13alpha),kmadd(J12L,kmul(J33L,PDstandardNth13alpha),kmadd(dJ123L,PDstandardNth1alpha,kmadd(J22L,kmul(J23L,PDstandardNth22alpha),kmadd(J23L,kmul(J32L,PDstandardNth23alpha),kmadd(J22L,kmul(J33L,PDstandardNth23alpha),kmadd(dJ223L,PDstandardNth2alpha,kmadd(J32L,kmul(J33L,PDstandardNth33alpha),kmul(dJ323L,PDstandardNth3alpha))))))))))));
      
      JacPDstandardNth23beta1 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11beta1),kmadd(J13L,kmul(J22L,PDstandardNth12beta1),kmadd(J12L,kmul(J23L,PDstandardNth12beta1),kmadd(J13L,kmul(J32L,PDstandardNth13beta1),kmadd(J12L,kmul(J33L,PDstandardNth13beta1),kmadd(dJ123L,PDstandardNth1beta1,kmadd(J22L,kmul(J23L,PDstandardNth22beta1),kmadd(J23L,kmul(J32L,PDstandardNth23beta1),kmadd(J22L,kmul(J33L,PDstandardNth23beta1),kmadd(dJ223L,PDstandardNth2beta1,kmadd(J32L,kmul(J33L,PDstandardNth33beta1),kmul(dJ323L,PDstandardNth3beta1))))))))))));
      
      JacPDstandardNth23beta2 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11beta2),kmadd(J13L,kmul(J22L,PDstandardNth12beta2),kmadd(J12L,kmul(J23L,PDstandardNth12beta2),kmadd(J13L,kmul(J32L,PDstandardNth13beta2),kmadd(J12L,kmul(J33L,PDstandardNth13beta2),kmadd(dJ123L,PDstandardNth1beta2,kmadd(J22L,kmul(J23L,PDstandardNth22beta2),kmadd(J23L,kmul(J32L,PDstandardNth23beta2),kmadd(J22L,kmul(J33L,PDstandardNth23beta2),kmadd(dJ223L,PDstandardNth2beta2,kmadd(J32L,kmul(J33L,PDstandardNth33beta2),kmul(dJ323L,PDstandardNth3beta2))))))))))));
      
      JacPDstandardNth23beta3 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11beta3),kmadd(J13L,kmul(J22L,PDstandardNth12beta3),kmadd(J12L,kmul(J23L,PDstandardNth12beta3),kmadd(J13L,kmul(J32L,PDstandardNth13beta3),kmadd(J12L,kmul(J33L,PDstandardNth13beta3),kmadd(dJ123L,PDstandardNth1beta3,kmadd(J22L,kmul(J23L,PDstandardNth22beta3),kmadd(J23L,kmul(J32L,PDstandardNth23beta3),kmadd(J22L,kmul(J33L,PDstandardNth23beta3),kmadd(dJ223L,PDstandardNth2beta3,kmadd(J32L,kmul(J33L,PDstandardNth33beta3),kmul(dJ323L,PDstandardNth3beta3))))))))))));
      
      JacPDstandardNth31alpha = 
        kmadd(J11L,kmul(J13L,PDstandardNth11alpha),kmadd(J13L,kmul(J21L,PDstandardNth12alpha),kmadd(J11L,kmul(J23L,PDstandardNth12alpha),kmadd(J13L,kmul(J31L,PDstandardNth13alpha),kmadd(J11L,kmul(J33L,PDstandardNth13alpha),kmadd(dJ113L,PDstandardNth1alpha,kmadd(J21L,kmul(J23L,PDstandardNth22alpha),kmadd(J23L,kmul(J31L,PDstandardNth23alpha),kmadd(J21L,kmul(J33L,PDstandardNth23alpha),kmadd(dJ213L,PDstandardNth2alpha,kmadd(J31L,kmul(J33L,PDstandardNth33alpha),kmul(dJ313L,PDstandardNth3alpha))))))))))));
      
      JacPDstandardNth31beta1 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11beta1),kmadd(J13L,kmul(J21L,PDstandardNth12beta1),kmadd(J11L,kmul(J23L,PDstandardNth12beta1),kmadd(J13L,kmul(J31L,PDstandardNth13beta1),kmadd(J11L,kmul(J33L,PDstandardNth13beta1),kmadd(dJ113L,PDstandardNth1beta1,kmadd(J21L,kmul(J23L,PDstandardNth22beta1),kmadd(J23L,kmul(J31L,PDstandardNth23beta1),kmadd(J21L,kmul(J33L,PDstandardNth23beta1),kmadd(dJ213L,PDstandardNth2beta1,kmadd(J31L,kmul(J33L,PDstandardNth33beta1),kmul(dJ313L,PDstandardNth3beta1))))))))))));
      
      JacPDstandardNth31beta2 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11beta2),kmadd(J13L,kmul(J21L,PDstandardNth12beta2),kmadd(J11L,kmul(J23L,PDstandardNth12beta2),kmadd(J13L,kmul(J31L,PDstandardNth13beta2),kmadd(J11L,kmul(J33L,PDstandardNth13beta2),kmadd(dJ113L,PDstandardNth1beta2,kmadd(J21L,kmul(J23L,PDstandardNth22beta2),kmadd(J23L,kmul(J31L,PDstandardNth23beta2),kmadd(J21L,kmul(J33L,PDstandardNth23beta2),kmadd(dJ213L,PDstandardNth2beta2,kmadd(J31L,kmul(J33L,PDstandardNth33beta2),kmul(dJ313L,PDstandardNth3beta2))))))))))));
      
      JacPDstandardNth31beta3 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11beta3),kmadd(J13L,kmul(J21L,PDstandardNth12beta3),kmadd(J11L,kmul(J23L,PDstandardNth12beta3),kmadd(J13L,kmul(J31L,PDstandardNth13beta3),kmadd(J11L,kmul(J33L,PDstandardNth13beta3),kmadd(dJ113L,PDstandardNth1beta3,kmadd(J21L,kmul(J23L,PDstandardNth22beta3),kmadd(J23L,kmul(J31L,PDstandardNth23beta3),kmadd(J21L,kmul(J33L,PDstandardNth23beta3),kmadd(dJ213L,PDstandardNth2beta3,kmadd(J31L,kmul(J33L,PDstandardNth33beta3),kmul(dJ313L,PDstandardNth3beta3))))))))))));
      
      JacPDstandardNth32alpha = 
        kmadd(J12L,kmul(J13L,PDstandardNth11alpha),kmadd(J13L,kmul(J22L,PDstandardNth12alpha),kmadd(J12L,kmul(J23L,PDstandardNth12alpha),kmadd(J13L,kmul(J32L,PDstandardNth13alpha),kmadd(J12L,kmul(J33L,PDstandardNth13alpha),kmadd(dJ123L,PDstandardNth1alpha,kmadd(J22L,kmul(J23L,PDstandardNth22alpha),kmadd(J23L,kmul(J32L,PDstandardNth23alpha),kmadd(J22L,kmul(J33L,PDstandardNth23alpha),kmadd(dJ223L,PDstandardNth2alpha,kmadd(J32L,kmul(J33L,PDstandardNth33alpha),kmul(dJ323L,PDstandardNth3alpha))))))))))));
      
      JacPDstandardNth32beta1 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11beta1),kmadd(J13L,kmul(J22L,PDstandardNth12beta1),kmadd(J12L,kmul(J23L,PDstandardNth12beta1),kmadd(J13L,kmul(J32L,PDstandardNth13beta1),kmadd(J12L,kmul(J33L,PDstandardNth13beta1),kmadd(dJ123L,PDstandardNth1beta1,kmadd(J22L,kmul(J23L,PDstandardNth22beta1),kmadd(J23L,kmul(J32L,PDstandardNth23beta1),kmadd(J22L,kmul(J33L,PDstandardNth23beta1),kmadd(dJ223L,PDstandardNth2beta1,kmadd(J32L,kmul(J33L,PDstandardNth33beta1),kmul(dJ323L,PDstandardNth3beta1))))))))))));
      
      JacPDstandardNth32beta2 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11beta2),kmadd(J13L,kmul(J22L,PDstandardNth12beta2),kmadd(J12L,kmul(J23L,PDstandardNth12beta2),kmadd(J13L,kmul(J32L,PDstandardNth13beta2),kmadd(J12L,kmul(J33L,PDstandardNth13beta2),kmadd(dJ123L,PDstandardNth1beta2,kmadd(J22L,kmul(J23L,PDstandardNth22beta2),kmadd(J23L,kmul(J32L,PDstandardNth23beta2),kmadd(J22L,kmul(J33L,PDstandardNth23beta2),kmadd(dJ223L,PDstandardNth2beta2,kmadd(J32L,kmul(J33L,PDstandardNth33beta2),kmul(dJ323L,PDstandardNth3beta2))))))))))));
      
      JacPDstandardNth32beta3 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11beta3),kmadd(J13L,kmul(J22L,PDstandardNth12beta3),kmadd(J12L,kmul(J23L,PDstandardNth12beta3),kmadd(J13L,kmul(J32L,PDstandardNth13beta3),kmadd(J12L,kmul(J33L,PDstandardNth13beta3),kmadd(dJ123L,PDstandardNth1beta3,kmadd(J22L,kmul(J23L,PDstandardNth22beta3),kmadd(J23L,kmul(J32L,PDstandardNth23beta3),kmadd(J22L,kmul(J33L,PDstandardNth23beta3),kmadd(dJ223L,PDstandardNth2beta3,kmadd(J32L,kmul(J33L,PDstandardNth33beta3),kmul(dJ323L,PDstandardNth3beta3))))))))))));
    }
    else
    {
      JacPDstandardNth1alpha = PDstandardNth1alpha;
      
      JacPDstandardNth1beta1 = PDstandardNth1beta1;
      
      JacPDstandardNth1beta2 = PDstandardNth1beta2;
      
      JacPDstandardNth1beta3 = PDstandardNth1beta3;
      
      JacPDstandardNth1gt11 = PDstandardNth1gt11;
      
      JacPDstandardNth1gt12 = PDstandardNth1gt12;
      
      JacPDstandardNth1gt13 = PDstandardNth1gt13;
      
      JacPDstandardNth1gt22 = PDstandardNth1gt22;
      
      JacPDstandardNth1gt23 = PDstandardNth1gt23;
      
      JacPDstandardNth1gt33 = PDstandardNth1gt33;
      
      JacPDstandardNth1phi = PDstandardNth1phi;
      
      JacPDstandardNth1trK = PDstandardNth1trK;
      
      JacPDstandardNth2alpha = PDstandardNth2alpha;
      
      JacPDstandardNth2beta1 = PDstandardNth2beta1;
      
      JacPDstandardNth2beta2 = PDstandardNth2beta2;
      
      JacPDstandardNth2beta3 = PDstandardNth2beta3;
      
      JacPDstandardNth2gt11 = PDstandardNth2gt11;
      
      JacPDstandardNth2gt12 = PDstandardNth2gt12;
      
      JacPDstandardNth2gt13 = PDstandardNth2gt13;
      
      JacPDstandardNth2gt22 = PDstandardNth2gt22;
      
      JacPDstandardNth2gt23 = PDstandardNth2gt23;
      
      JacPDstandardNth2gt33 = PDstandardNth2gt33;
      
      JacPDstandardNth2phi = PDstandardNth2phi;
      
      JacPDstandardNth2trK = PDstandardNth2trK;
      
      JacPDstandardNth3alpha = PDstandardNth3alpha;
      
      JacPDstandardNth3beta1 = PDstandardNth3beta1;
      
      JacPDstandardNth3beta2 = PDstandardNth3beta2;
      
      JacPDstandardNth3beta3 = PDstandardNth3beta3;
      
      JacPDstandardNth3gt11 = PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = PDstandardNth3gt33;
      
      JacPDstandardNth3phi = PDstandardNth3phi;
      
      JacPDstandardNth3trK = PDstandardNth3trK;
      
      JacPDstandardNth11alpha = PDstandardNth11alpha;
      
      JacPDstandardNth11beta1 = PDstandardNth11beta1;
      
      JacPDstandardNth11beta2 = PDstandardNth11beta2;
      
      JacPDstandardNth11beta3 = PDstandardNth11beta3;
      
      JacPDstandardNth22alpha = PDstandardNth22alpha;
      
      JacPDstandardNth22beta1 = PDstandardNth22beta1;
      
      JacPDstandardNth22beta2 = PDstandardNth22beta2;
      
      JacPDstandardNth22beta3 = PDstandardNth22beta3;
      
      JacPDstandardNth33alpha = PDstandardNth33alpha;
      
      JacPDstandardNth33beta1 = PDstandardNth33beta1;
      
      JacPDstandardNth33beta2 = PDstandardNth33beta2;
      
      JacPDstandardNth33beta3 = PDstandardNth33beta3;
      
      JacPDstandardNth12alpha = PDstandardNth12alpha;
      
      JacPDstandardNth12beta1 = PDstandardNth12beta1;
      
      JacPDstandardNth12beta2 = PDstandardNth12beta2;
      
      JacPDstandardNth12beta3 = PDstandardNth12beta3;
      
      JacPDstandardNth13alpha = PDstandardNth13alpha;
      
      JacPDstandardNth13beta1 = PDstandardNth13beta1;
      
      JacPDstandardNth13beta2 = PDstandardNth13beta2;
      
      JacPDstandardNth13beta3 = PDstandardNth13beta3;
      
      JacPDstandardNth21alpha = PDstandardNth12alpha;
      
      JacPDstandardNth21beta1 = PDstandardNth12beta1;
      
      JacPDstandardNth21beta2 = PDstandardNth12beta2;
      
      JacPDstandardNth21beta3 = PDstandardNth12beta3;
      
      JacPDstandardNth23alpha = PDstandardNth23alpha;
      
      JacPDstandardNth23beta1 = PDstandardNth23beta1;
      
      JacPDstandardNth23beta2 = PDstandardNth23beta2;
      
      JacPDstandardNth23beta3 = PDstandardNth23beta3;
      
      JacPDstandardNth31alpha = PDstandardNth13alpha;
      
      JacPDstandardNth31beta1 = PDstandardNth13beta1;
      
      JacPDstandardNth31beta2 = PDstandardNth13beta2;
      
      JacPDstandardNth31beta3 = PDstandardNth13beta3;
      
      JacPDstandardNth32alpha = PDstandardNth23alpha;
      
      JacPDstandardNth32beta1 = PDstandardNth23beta1;
      
      JacPDstandardNth32beta2 = PDstandardNth23beta2;
      
      JacPDstandardNth32beta3 = PDstandardNth23beta3;
    }
    
    CCTK_REAL_VEC detgt CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC gtu11 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt22L,gt33L,kmul(gt23L,gt23L)),detgt);
    
    CCTK_REAL_VEC gtu12 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt13L,gt23L,kmul(gt12L,gt33L)),detgt);
    
    CCTK_REAL_VEC gtu13 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt12L,gt23L,kmul(gt13L,gt22L)),detgt);
    
    CCTK_REAL_VEC gtu22 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt11L,gt33L,kmul(gt13L,gt13L)),detgt);
    
    CCTK_REAL_VEC gtu23 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt12L,gt13L,kmul(gt11L,gt23L)),detgt);
    
    CCTK_REAL_VEC gtu33 CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kmsub(gt11L,gt22L,kmul(gt12L,gt12L)),detgt);
    
    CCTK_REAL_VEC Gtl111 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth1gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl112 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth2gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl113 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth3gt11,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl122 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmsub(JacPDstandardNth2gt12,ToReal(2),JacPDstandardNth1gt22));
    
    CCTK_REAL_VEC Gtl123 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kadd(JacPDstandardNth2gt13,ksub(JacPDstandardNth3gt12,JacPDstandardNth1gt23)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl133 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmsub(JacPDstandardNth3gt13,ToReal(2),JacPDstandardNth1gt33));
    
    CCTK_REAL_VEC Gtl211 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmsub(JacPDstandardNth1gt12,ToReal(2),JacPDstandardNth2gt11));
    
    CCTK_REAL_VEC Gtl212 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth1gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl213 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth3gt12,JacPDstandardNth2gt13)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl222 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth2gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl223 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth3gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl233 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmsub(JacPDstandardNth3gt23,ToReal(2),JacPDstandardNth2gt33));
    
    CCTK_REAL_VEC Gtl311 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmsub(JacPDstandardNth1gt13,ToReal(2),JacPDstandardNth3gt11));
    
    CCTK_REAL_VEC Gtl312 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth2gt13,JacPDstandardNth3gt12)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl313 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth1gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl322 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ToReal(0.5),kmsub(JacPDstandardNth2gt23,ToReal(2),JacPDstandardNth3gt22));
    
    CCTK_REAL_VEC Gtl323 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth2gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl333 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth3gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gt111 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu11,kmadd(Gtl211,gtu12,kmul(Gtl311,gtu13)));
    
    CCTK_REAL_VEC Gt211 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu12,kmadd(Gtl211,gtu22,kmul(Gtl311,gtu23)));
    
    CCTK_REAL_VEC Gt311 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu13,kmadd(Gtl211,gtu23,kmul(Gtl311,gtu33)));
    
    CCTK_REAL_VEC Gt112 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu11,kmadd(Gtl212,gtu12,kmul(Gtl312,gtu13)));
    
    CCTK_REAL_VEC Gt212 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu12,kmadd(Gtl212,gtu22,kmul(Gtl312,gtu23)));
    
    CCTK_REAL_VEC Gt312 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu13,kmadd(Gtl212,gtu23,kmul(Gtl312,gtu33)));
    
    CCTK_REAL_VEC Gt113 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu11,kmadd(Gtl213,gtu12,kmul(Gtl313,gtu13)));
    
    CCTK_REAL_VEC Gt213 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu12,kmadd(Gtl213,gtu22,kmul(Gtl313,gtu23)));
    
    CCTK_REAL_VEC Gt313 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu13,kmadd(Gtl213,gtu23,kmul(Gtl313,gtu33)));
    
    CCTK_REAL_VEC Gt122 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl122,gtu11,kmadd(Gtl222,gtu12,kmul(Gtl322,gtu13)));
    
    CCTK_REAL_VEC Gt222 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl122,gtu12,kmadd(Gtl222,gtu22,kmul(Gtl322,gtu23)));
    
    CCTK_REAL_VEC Gt322 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl122,gtu13,kmadd(Gtl222,gtu23,kmul(Gtl322,gtu33)));
    
    CCTK_REAL_VEC Gt123 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl123,gtu11,kmadd(Gtl223,gtu12,kmul(Gtl323,gtu13)));
    
    CCTK_REAL_VEC Gt223 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl123,gtu12,kmadd(Gtl223,gtu22,kmul(Gtl323,gtu23)));
    
    CCTK_REAL_VEC Gt323 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl123,gtu13,kmadd(Gtl223,gtu23,kmul(Gtl323,gtu33)));
    
    CCTK_REAL_VEC Gt133 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl133,gtu11,kmadd(Gtl233,gtu12,kmul(Gtl333,gtu13)));
    
    CCTK_REAL_VEC Gt233 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl133,gtu12,kmadd(Gtl233,gtu22,kmul(Gtl333,gtu23)));
    
    CCTK_REAL_VEC Gt333 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl133,gtu13,kmadd(Gtl233,gtu23,kmul(Gtl333,gtu33)));
    
    CCTK_REAL_VEC Xtn1 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt111,gtu11,kmadd(Gt122,gtu22,kmadd(Gt133,gtu33,kmadd(Gt112,kmul(gtu12,ToReal(2)),kmadd(Gt113,kmul(gtu13,ToReal(2)),kmul(Gt123,kmul(gtu23,ToReal(2))))))));
    
    CCTK_REAL_VEC Xtn2 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt211,gtu11,kmadd(Gt222,gtu22,kmadd(Gt233,gtu33,kmadd(Gt212,kmul(gtu12,ToReal(2)),kmadd(Gt213,kmul(gtu13,ToReal(2)),kmul(Gt223,kmul(gtu23,ToReal(2))))))));
    
    CCTK_REAL_VEC Xtn3 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt311,gtu11,kmadd(Gt322,gtu22,kmadd(Gt333,gtu33,kmadd(Gt312,kmul(gtu12,ToReal(2)),kmadd(Gt313,kmul(gtu13,ToReal(2)),kmul(Gt323,kmul(gtu23,ToReal(2))))))));
    
    CCTK_REAL_VEC e4phi CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod == 
      1,kdiv(ToReal(1),kmul(phiL,phiL)),kexp(kmul(phiL,ToReal(4))));
    
    CCTK_REAL_VEC em4phi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),e4phi);
    
    CCTK_REAL_VEC fac1 CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod == 
      1,kdiv(ToReal(-0.5),phiL),ToReal(1));
    
    CCTK_REAL_VEC cdphi1 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,JacPDstandardNth1phi);
    
    CCTK_REAL_VEC cdphi2 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,JacPDstandardNth2phi);
    
    CCTK_REAL_VEC cdphi3 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,JacPDstandardNth3phi);
    
    CCTK_REAL_VEC Atm11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At11L,gtu11,kmadd(At12L,gtu12,kmul(At13L,gtu13)));
    
    CCTK_REAL_VEC Atm21 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At11L,gtu12,kmadd(At12L,gtu22,kmul(At13L,gtu23)));
    
    CCTK_REAL_VEC Atm31 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At11L,gtu13,kmadd(At12L,gtu23,kmul(At13L,gtu33)));
    
    CCTK_REAL_VEC Atm12 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At12L,gtu11,kmadd(At22L,gtu12,kmul(At23L,gtu13)));
    
    CCTK_REAL_VEC Atm22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At12L,gtu12,kmadd(At22L,gtu22,kmul(At23L,gtu23)));
    
    CCTK_REAL_VEC Atm32 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At12L,gtu13,kmadd(At22L,gtu23,kmul(At23L,gtu33)));
    
    CCTK_REAL_VEC Atm13 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At13L,gtu11,kmadd(At23L,gtu12,kmul(At33L,gtu13)));
    
    CCTK_REAL_VEC Atm23 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At13L,gtu12,kmadd(At23L,gtu22,kmul(At33L,gtu23)));
    
    CCTK_REAL_VEC Atm33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(At13L,gtu13,kmadd(At23L,gtu23,kmul(At33L,gtu33)));
    
    CCTK_REAL_VEC Atu11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm11,gtu11,kmadd(Atm12,gtu12,kmul(Atm13,gtu13)));
    
    CCTK_REAL_VEC Atu12 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm11,gtu12,kmadd(Atm12,gtu22,kmul(Atm13,gtu23)));
    
    CCTK_REAL_VEC Atu13 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm11,gtu13,kmadd(Atm12,gtu23,kmul(Atm13,gtu33)));
    
    CCTK_REAL_VEC Atu22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm21,gtu12,kmadd(Atm22,gtu22,kmul(Atm23,gtu23)));
    
    CCTK_REAL_VEC Atu23 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm21,gtu13,kmadd(Atm22,gtu23,kmul(Atm23,gtu33)));
    
    CCTK_REAL_VEC Atu33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Atm31,gtu13,kmadd(Atm32,gtu23,kmul(Atm33,gtu33)));
    
    CCTK_REAL_VEC rho CCTK_ATTRIBUTE_UNUSED = 
      kdiv(kadd(eTttL,kmadd(beta1L,kmadd(beta1L,eTxxL,kmadd(beta2L,eTxyL,kmul(beta3L,eTxzL))),kmadd(beta2L,kmadd(beta1L,eTxyL,kmadd(beta2L,eTyyL,kmul(beta3L,eTyzL))),kmadd(beta3L,kmadd(beta1L,eTxzL,kmadd(beta2L,eTyzL,kmul(beta3L,eTzzL))),kmul(kmadd(beta1L,eTtxL,kmadd(beta2L,eTtyL,kmul(beta3L,eTtzL))),ToReal(-2)))))),kmul(alphaL,alphaL));
    
    CCTK_REAL_VEC S1 CCTK_ATTRIBUTE_UNUSED = 
      kneg(kdiv(ksub(eTtxL,kmadd(beta1L,eTxxL,kmadd(beta3L,eTxzL,kmul(beta2L,eTxyL)))),alphaL));
    
    CCTK_REAL_VEC S2 CCTK_ATTRIBUTE_UNUSED = 
      kneg(kdiv(ksub(eTtyL,kmadd(beta1L,eTxyL,kmadd(beta3L,eTyzL,kmul(beta2L,eTyyL)))),alphaL));
    
    CCTK_REAL_VEC S3 CCTK_ATTRIBUTE_UNUSED = 
      kneg(kdiv(ksub(eTtzL,kmadd(beta1L,eTxzL,kmadd(beta3L,eTzzL,kmul(beta2L,eTyzL)))),alphaL));
    
    CCTK_REAL_VEC trS CCTK_ATTRIBUTE_UNUSED = 
      kmul(em4phi,kmadd(eTxxL,gtu11,kmadd(eTyyL,gtu22,kmadd(eTzzL,gtu33,kmadd(eTxyL,kmul(gtu12,ToReal(2)),kmadd(eTxzL,kmul(gtu13,ToReal(2)),kmul(eTyzL,kmul(gtu23,ToReal(2)))))))));
    
    CCTK_REAL_VEC phirhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(IfThen(conformalMethod == 
      1,kmul(phiL,ToReal(0.333333333333333333333333333333)),ToReal(-0.166666666666666666666666666667)),kmsub(alphaL,trKL,kadd(kadd(JacPDstandardNth3beta3,JacPDstandardNth2beta2),JacPDstandardNth1beta1)));
    
    CCTK_REAL_VEC gt11rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmul(At11L,ToReal(-2)),kmadd(gt11L,kmul(kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667)),kmadd(gt11L,kmul(JacPDstandardNth1beta1,ToReal(2)),kmadd(gt12L,kmul(JacPDstandardNth1beta2,ToReal(2)),kmul(gt13L,kmul(JacPDstandardNth1beta3,ToReal(2)))))));
    
    CCTK_REAL_VEC gt12rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gt12L,JacPDstandardNth1beta1,kmadd(gt22L,JacPDstandardNth1beta2,kmadd(gt23L,JacPDstandardNth1beta3,kmadd(gt11L,JacPDstandardNth2beta1,kmadd(gt12L,JacPDstandardNth2beta2,kmadd(gt13L,JacPDstandardNth2beta3,kmadd(alphaL,kmul(At12L,ToReal(-2)),kmul(gt12L,kmul(kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667))))))))));
    
    CCTK_REAL_VEC gt13rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gt13L,JacPDstandardNth1beta1,kmadd(gt23L,JacPDstandardNth1beta2,kmadd(gt33L,JacPDstandardNth1beta3,kmadd(gt11L,JacPDstandardNth3beta1,kmadd(gt12L,JacPDstandardNth3beta2,kmadd(gt13L,JacPDstandardNth3beta3,kmadd(alphaL,kmul(At13L,ToReal(-2)),kmul(gt13L,kmul(kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667))))))))));
    
    CCTK_REAL_VEC gt22rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmul(At22L,ToReal(-2)),kmadd(gt22L,kmul(kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667)),kmadd(gt12L,kmul(JacPDstandardNth2beta1,ToReal(2)),kmadd(gt22L,kmul(JacPDstandardNth2beta2,ToReal(2)),kmul(gt23L,kmul(JacPDstandardNth2beta3,ToReal(2)))))));
    
    CCTK_REAL_VEC gt23rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gt13L,JacPDstandardNth2beta1,kmadd(gt23L,JacPDstandardNth2beta2,kmadd(gt33L,JacPDstandardNth2beta3,kmadd(gt12L,JacPDstandardNth3beta1,kmadd(gt22L,JacPDstandardNth3beta2,kmadd(gt23L,JacPDstandardNth3beta3,kmadd(alphaL,kmul(At23L,ToReal(-2)),kmul(gt23L,kmul(kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667))))))))));
    
    CCTK_REAL_VEC gt33rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmul(At33L,ToReal(-2)),kmadd(gt33L,kmul(kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)),ToReal(-0.666666666666666666666666666667)),kmadd(gt13L,kmul(JacPDstandardNth3beta1,ToReal(2)),kmadd(gt23L,kmul(JacPDstandardNth3beta2,ToReal(2)),kmul(gt33L,kmul(JacPDstandardNth3beta3,ToReal(2)))))));
    
    CCTK_REAL_VEC dotXt1 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gtu11,JacPDstandardNth11beta1,kmadd(gtu12,JacPDstandardNth12beta1,kmadd(gtu13,JacPDstandardNth13beta1,kmadd(gtu12,JacPDstandardNth21beta1,kmadd(gtu22,JacPDstandardNth22beta1,kmadd(gtu23,JacPDstandardNth23beta1,kmadd(gtu13,JacPDstandardNth31beta1,kmadd(gtu23,JacPDstandardNth32beta1,kmadd(gtu33,JacPDstandardNth33beta1,kmadd(kmadd(Atu11,JacPDstandardNth1alpha,kmadd(Atu12,JacPDstandardNth2alpha,kmul(Atu13,JacPDstandardNth3alpha))),ToReal(-2),knmsub(JacPDstandardNth1beta1,Xtn1,knmsub(JacPDstandardNth2beta1,Xtn2,knmsub(JacPDstandardNth3beta1,Xtn3,kmadd(kmadd(gtu11,kadd(JacPDstandardNth11beta1,kadd(JacPDstandardNth12beta2,JacPDstandardNth13beta3)),kmadd(gtu12,kadd(JacPDstandardNth21beta1,kadd(JacPDstandardNth22beta2,JacPDstandardNth23beta3)),kmul(gtu13,kadd(JacPDstandardNth31beta1,kadd(JacPDstandardNth32beta2,JacPDstandardNth33beta3))))),ToReal(0.333333333333333333333333333333),kmadd(kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)),kmul(Xtn1,ToReal(0.666666666666666666666666666667)),kmadd(alphaL,kmul(ToReal(2),kmadd(Atu11,Gt111,kmadd(Atu22,Gt122,kmadd(Atu33,Gt133,kmadd(kmadd(gtu11,JacPDstandardNth1trK,kmadd(gtu12,JacPDstandardNth2trK,kmul(gtu13,JacPDstandardNth3trK))),ToReal(-0.666666666666666666666666666667),kmadd(Atu12,kmul(Gt112,ToReal(2)),kmadd(Atu13,kmul(Gt113,ToReal(2)),kmadd(Atu23,kmul(Gt123,ToReal(2)),kmul(kmadd(Atu11,cdphi1,kmadd(Atu12,cdphi2,kmul(Atu13,cdphi3))),ToReal(6)))))))))),kmul(alphaL,kmul(kmadd(gtu11,S1,kmadd(gtu12,S2,kmul(gtu13,S3))),ToReal(-50.2654824574366918154022941325)))))))))))))))))));
    
    CCTK_REAL_VEC dotXt2 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gtu11,JacPDstandardNth11beta2,kmadd(gtu12,JacPDstandardNth12beta2,kmadd(gtu13,JacPDstandardNth13beta2,kmadd(gtu12,JacPDstandardNth21beta2,kmadd(gtu22,JacPDstandardNth22beta2,kmadd(gtu23,JacPDstandardNth23beta2,kmadd(gtu13,JacPDstandardNth31beta2,kmadd(gtu23,JacPDstandardNth32beta2,kmadd(gtu33,JacPDstandardNth33beta2,kmadd(kmadd(Atu12,JacPDstandardNth1alpha,kmadd(Atu22,JacPDstandardNth2alpha,kmul(Atu23,JacPDstandardNth3alpha))),ToReal(-2),knmsub(JacPDstandardNth1beta2,Xtn1,knmsub(JacPDstandardNth2beta2,Xtn2,knmsub(JacPDstandardNth3beta2,Xtn3,kmadd(kmadd(gtu12,kadd(JacPDstandardNth11beta1,kadd(JacPDstandardNth12beta2,JacPDstandardNth13beta3)),kmadd(gtu22,kadd(JacPDstandardNth21beta1,kadd(JacPDstandardNth22beta2,JacPDstandardNth23beta3)),kmul(gtu23,kadd(JacPDstandardNth31beta1,kadd(JacPDstandardNth32beta2,JacPDstandardNth33beta3))))),ToReal(0.333333333333333333333333333333),kmadd(kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)),kmul(Xtn2,ToReal(0.666666666666666666666666666667)),kmadd(alphaL,kmul(ToReal(2),kmadd(Atu11,Gt211,kmadd(Atu22,Gt222,kmadd(Atu33,Gt233,kmadd(kmadd(gtu12,JacPDstandardNth1trK,kmadd(gtu22,JacPDstandardNth2trK,kmul(gtu23,JacPDstandardNth3trK))),ToReal(-0.666666666666666666666666666667),kmadd(Atu12,kmul(Gt212,ToReal(2)),kmadd(Atu13,kmul(Gt213,ToReal(2)),kmadd(Atu23,kmul(Gt223,ToReal(2)),kmul(kmadd(Atu12,cdphi1,kmadd(Atu22,cdphi2,kmul(Atu23,cdphi3))),ToReal(6)))))))))),kmul(alphaL,kmul(kmadd(gtu12,S1,kmadd(gtu22,S2,kmul(gtu23,S3))),ToReal(-50.2654824574366918154022941325)))))))))))))))))));
    
    CCTK_REAL_VEC dotXt3 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gtu11,JacPDstandardNth11beta3,kmadd(gtu12,JacPDstandardNth12beta3,kmadd(gtu13,JacPDstandardNth13beta3,kmadd(gtu12,JacPDstandardNth21beta3,kmadd(gtu22,JacPDstandardNth22beta3,kmadd(gtu23,JacPDstandardNth23beta3,kmadd(gtu13,JacPDstandardNth31beta3,kmadd(gtu23,JacPDstandardNth32beta3,kmadd(gtu33,JacPDstandardNth33beta3,kmadd(kmadd(Atu13,JacPDstandardNth1alpha,kmadd(Atu23,JacPDstandardNth2alpha,kmul(Atu33,JacPDstandardNth3alpha))),ToReal(-2),knmsub(JacPDstandardNth1beta3,Xtn1,knmsub(JacPDstandardNth2beta3,Xtn2,knmsub(JacPDstandardNth3beta3,Xtn3,kmadd(kmadd(gtu13,kadd(JacPDstandardNth11beta1,kadd(JacPDstandardNth12beta2,JacPDstandardNth13beta3)),kmadd(gtu23,kadd(JacPDstandardNth21beta1,kadd(JacPDstandardNth22beta2,JacPDstandardNth23beta3)),kmul(gtu33,kadd(JacPDstandardNth31beta1,kadd(JacPDstandardNth32beta2,JacPDstandardNth33beta3))))),ToReal(0.333333333333333333333333333333),kmadd(kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3)),kmul(Xtn3,ToReal(0.666666666666666666666666666667)),kmadd(alphaL,kmul(ToReal(2),kmadd(Atu11,Gt311,kmadd(Atu22,Gt322,kmadd(Atu33,Gt333,kmadd(kmadd(gtu13,JacPDstandardNth1trK,kmadd(gtu23,JacPDstandardNth2trK,kmul(gtu33,JacPDstandardNth3trK))),ToReal(-0.666666666666666666666666666667),kmadd(Atu12,kmul(Gt312,ToReal(2)),kmadd(Atu13,kmul(Gt313,ToReal(2)),kmadd(Atu23,kmul(Gt323,ToReal(2)),kmul(kmadd(Atu13,cdphi1,kmadd(Atu23,cdphi2,kmul(Atu33,cdphi3))),ToReal(6)))))))))),kmul(alphaL,kmul(kmadd(gtu13,S1,kmadd(gtu23,S2,kmul(gtu33,S3))),ToReal(-50.2654824574366918154022941325)))))))))))))))))));
    
    CCTK_REAL_VEC Xt1rhsL CCTK_ATTRIBUTE_UNUSED = dotXt1;
    
    CCTK_REAL_VEC Xt2rhsL CCTK_ATTRIBUTE_UNUSED = dotXt2;
    
    CCTK_REAL_VEC Xt3rhsL CCTK_ATTRIBUTE_UNUSED = dotXt3;
    
    CCTK_REAL_VEC dottrK CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmadd(Atm11,Atm11,kmadd(Atm22,Atm22,kmadd(Atm33,Atm33,kmadd(kmul(trKL,trKL),ToReal(0.333333333333333333333333333333),kmadd(Atm12,kmul(Atm21,ToReal(2)),kmadd(Atm13,kmul(Atm31,ToReal(2)),kmul(Atm23,kmul(Atm32,ToReal(2))))))))),kmsub(alphaL,kmul(kadd(rho,trS),ToReal(12.5663706143591729538505735331)),kmul(em4phi,knmsub(JacPDstandardNth1alpha,Xtn1,knmsub(JacPDstandardNth2alpha,Xtn2,knmsub(JacPDstandardNth3alpha,Xtn3,kmadd(gtu11,kmadd(cdphi1,kmul(JacPDstandardNth1alpha,ToReal(2)),JacPDstandardNth11alpha),kmadd(gtu12,kmadd(cdphi2,kmul(JacPDstandardNth1alpha,ToReal(2)),JacPDstandardNth21alpha),kmadd(gtu13,kmadd(cdphi3,kmul(JacPDstandardNth1alpha,ToReal(2)),JacPDstandardNth31alpha),kmadd(gtu12,kmadd(cdphi1,kmul(JacPDstandardNth2alpha,ToReal(2)),JacPDstandardNth12alpha),kmadd(gtu22,kmadd(cdphi2,kmul(JacPDstandardNth2alpha,ToReal(2)),JacPDstandardNth22alpha),kmadd(gtu23,kmadd(cdphi3,kmul(JacPDstandardNth2alpha,ToReal(2)),JacPDstandardNth32alpha),kmadd(gtu13,kmadd(cdphi1,kmul(JacPDstandardNth3alpha,ToReal(2)),JacPDstandardNth13alpha),kmadd(gtu23,kmadd(cdphi2,kmul(JacPDstandardNth3alpha,ToReal(2)),JacPDstandardNth23alpha),kmul(gtu33,kmadd(cdphi3,kmul(JacPDstandardNth3alpha,ToReal(2)),JacPDstandardNth33alpha))))))))))))))));
    
    CCTK_REAL_VEC trKrhsL CCTK_ATTRIBUTE_UNUSED = dottrK;
    
    CCTK_REAL_VEC alpharhsL CCTK_ATTRIBUTE_UNUSED = 
      kneg(kmul(kpow(alphaL,harmonicN),kmul(ToReal(harmonicF),kmadd(AL,ToReal(LapseACoeff),kmul(kmadd(kadd(alphaL,ToReal(-1)),ToReal(AlphaDriver),trKL),ToReal(1 
      - LapseACoeff))))));
    
    CCTK_REAL_VEC ArhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(knmsub(AL,ToReal(AlphaDriver),dottrK),ToReal(LapseACoeff));
    
    CCTK_REAL_VEC eta CCTK_ATTRIBUTE_UNUSED = 
      kdiv(ToReal(SpatialBetaDriverRadius),kfmax(rL,ToReal(SpatialBetaDriverRadius)));
    
    CCTK_REAL_VEC theta CCTK_ATTRIBUTE_UNUSED = 
      kfmin(ToReal(1),kexp(knmsub(rL,ToReal(ScalarINV(SpatialShiftGammaCoeffRadius)),ToReal(1))));
    
    CCTK_REAL_VEC Ddetgt1 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gtu11,JacPDstandardNth1gt11,kmadd(gtu22,JacPDstandardNth1gt22,kmadd(gtu33,JacPDstandardNth1gt33,kmadd(gtu12,kmul(JacPDstandardNth1gt12,ToReal(2)),kmadd(gtu13,kmul(JacPDstandardNth1gt13,ToReal(2)),kmul(gtu23,kmul(JacPDstandardNth1gt23,ToReal(2))))))));
    
    CCTK_REAL_VEC Ddetgt2 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gtu11,JacPDstandardNth2gt11,kmadd(gtu22,JacPDstandardNth2gt22,kmadd(gtu33,JacPDstandardNth2gt33,kmadd(gtu12,kmul(JacPDstandardNth2gt12,ToReal(2)),kmadd(gtu13,kmul(JacPDstandardNth2gt13,ToReal(2)),kmul(gtu23,kmul(JacPDstandardNth2gt23,ToReal(2))))))));
    
    CCTK_REAL_VEC Ddetgt3 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(gtu11,JacPDstandardNth3gt11,kmadd(gtu22,JacPDstandardNth3gt22,kmadd(gtu33,JacPDstandardNth3gt33,kmadd(gtu12,kmul(JacPDstandardNth3gt12,ToReal(2)),kmadd(gtu13,kmul(JacPDstandardNth3gt13,ToReal(2)),kmul(gtu23,kmul(JacPDstandardNth3gt23,ToReal(2))))))));
    
    CCTK_REAL_VEC beta1rhsL CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC beta2rhsL CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC beta3rhsL CCTK_ATTRIBUTE_UNUSED;
    
    if (harmonicShift)
    {
      beta1rhsL = 
        kmul(alphaL,kmul(em4phi,kmul(ToReal(-0.5),kmadd(gtu11,kmadd(alphaL,kmul(JacPDstandardNth1phi,kmul(IfThen(conformalMethod 
        == 
        1,kdiv(ToReal(1),phiL),ToReal(-2)),ToReal(-2))),kmadd(alphaL,kmadd(kmadd(gtu11,JacPDstandardNth1gt11,kmadd(gtu12,JacPDstandardNth1gt12,kmadd(gtu13,JacPDstandardNth1gt13,kmadd(gtu12,JacPDstandardNth2gt11,kmadd(gtu22,JacPDstandardNth2gt12,kmadd(gtu23,JacPDstandardNth2gt13,kmadd(gtu13,JacPDstandardNth3gt11,kmadd(gtu23,JacPDstandardNth3gt12,kmul(gtu33,JacPDstandardNth3gt13))))))))),ToReal(-2),Ddetgt1),kmul(JacPDstandardNth1alpha,ToReal(2)))),kmadd(gtu12,kmadd(alphaL,kmul(JacPDstandardNth2phi,kmul(IfThen(conformalMethod 
        == 
        1,kdiv(ToReal(1),phiL),ToReal(-2)),ToReal(-2))),kmadd(alphaL,kmadd(kmadd(gtu11,JacPDstandardNth1gt12,kmadd(gtu12,JacPDstandardNth1gt22,kmadd(gtu13,JacPDstandardNth1gt23,kmadd(gtu12,JacPDstandardNth2gt12,kmadd(gtu22,JacPDstandardNth2gt22,kmadd(gtu23,JacPDstandardNth2gt23,kmadd(gtu13,JacPDstandardNth3gt12,kmadd(gtu23,JacPDstandardNth3gt22,kmul(gtu33,JacPDstandardNth3gt23))))))))),ToReal(-2),Ddetgt2),kmul(JacPDstandardNth2alpha,ToReal(2)))),kmul(gtu13,kmadd(alphaL,kmul(JacPDstandardNth3phi,kmul(IfThen(conformalMethod 
        == 
        1,kdiv(ToReal(1),phiL),ToReal(-2)),ToReal(-2))),kmadd(alphaL,kmadd(kmadd(gtu11,JacPDstandardNth1gt13,kmadd(gtu12,JacPDstandardNth1gt23,kmadd(gtu13,JacPDstandardNth1gt33,kmadd(gtu12,JacPDstandardNth2gt13,kmadd(gtu22,JacPDstandardNth2gt23,kmadd(gtu23,JacPDstandardNth2gt33,kmadd(gtu13,JacPDstandardNth3gt13,kmadd(gtu23,JacPDstandardNth3gt23,kmul(gtu33,JacPDstandardNth3gt33))))))))),ToReal(-2),Ddetgt3),kmul(JacPDstandardNth3alpha,ToReal(2))))))))));
      
      beta2rhsL = 
        kmul(alphaL,kmul(em4phi,kmul(ToReal(-0.5),kmadd(gtu12,kmadd(alphaL,kmul(JacPDstandardNth1phi,kmul(IfThen(conformalMethod 
        == 
        1,kdiv(ToReal(1),phiL),ToReal(-2)),ToReal(-2))),kmadd(alphaL,kmadd(kmadd(gtu11,JacPDstandardNth1gt11,kmadd(gtu12,JacPDstandardNth1gt12,kmadd(gtu13,JacPDstandardNth1gt13,kmadd(gtu12,JacPDstandardNth2gt11,kmadd(gtu22,JacPDstandardNth2gt12,kmadd(gtu23,JacPDstandardNth2gt13,kmadd(gtu13,JacPDstandardNth3gt11,kmadd(gtu23,JacPDstandardNth3gt12,kmul(gtu33,JacPDstandardNth3gt13))))))))),ToReal(-2),Ddetgt1),kmul(JacPDstandardNth1alpha,ToReal(2)))),kmadd(gtu22,kmadd(alphaL,kmul(JacPDstandardNth2phi,kmul(IfThen(conformalMethod 
        == 
        1,kdiv(ToReal(1),phiL),ToReal(-2)),ToReal(-2))),kmadd(alphaL,kmadd(kmadd(gtu11,JacPDstandardNth1gt12,kmadd(gtu12,JacPDstandardNth1gt22,kmadd(gtu13,JacPDstandardNth1gt23,kmadd(gtu12,JacPDstandardNth2gt12,kmadd(gtu22,JacPDstandardNth2gt22,kmadd(gtu23,JacPDstandardNth2gt23,kmadd(gtu13,JacPDstandardNth3gt12,kmadd(gtu23,JacPDstandardNth3gt22,kmul(gtu33,JacPDstandardNth3gt23))))))))),ToReal(-2),Ddetgt2),kmul(JacPDstandardNth2alpha,ToReal(2)))),kmul(gtu23,kmadd(alphaL,kmul(JacPDstandardNth3phi,kmul(IfThen(conformalMethod 
        == 
        1,kdiv(ToReal(1),phiL),ToReal(-2)),ToReal(-2))),kmadd(alphaL,kmadd(kmadd(gtu11,JacPDstandardNth1gt13,kmadd(gtu12,JacPDstandardNth1gt23,kmadd(gtu13,JacPDstandardNth1gt33,kmadd(gtu12,JacPDstandardNth2gt13,kmadd(gtu22,JacPDstandardNth2gt23,kmadd(gtu23,JacPDstandardNth2gt33,kmadd(gtu13,JacPDstandardNth3gt13,kmadd(gtu23,JacPDstandardNth3gt23,kmul(gtu33,JacPDstandardNth3gt33))))))))),ToReal(-2),Ddetgt3),kmul(JacPDstandardNth3alpha,ToReal(2))))))))));
      
      beta3rhsL = 
        kmul(alphaL,kmul(em4phi,kmul(ToReal(-0.5),kmadd(gtu13,kmadd(alphaL,kmul(JacPDstandardNth1phi,kmul(IfThen(conformalMethod 
        == 
        1,kdiv(ToReal(1),phiL),ToReal(-2)),ToReal(-2))),kmadd(alphaL,kmadd(kmadd(gtu11,JacPDstandardNth1gt11,kmadd(gtu12,JacPDstandardNth1gt12,kmadd(gtu13,JacPDstandardNth1gt13,kmadd(gtu12,JacPDstandardNth2gt11,kmadd(gtu22,JacPDstandardNth2gt12,kmadd(gtu23,JacPDstandardNth2gt13,kmadd(gtu13,JacPDstandardNth3gt11,kmadd(gtu23,JacPDstandardNth3gt12,kmul(gtu33,JacPDstandardNth3gt13))))))))),ToReal(-2),Ddetgt1),kmul(JacPDstandardNth1alpha,ToReal(2)))),kmadd(gtu23,kmadd(alphaL,kmul(JacPDstandardNth2phi,kmul(IfThen(conformalMethod 
        == 
        1,kdiv(ToReal(1),phiL),ToReal(-2)),ToReal(-2))),kmadd(alphaL,kmadd(kmadd(gtu11,JacPDstandardNth1gt12,kmadd(gtu12,JacPDstandardNth1gt22,kmadd(gtu13,JacPDstandardNth1gt23,kmadd(gtu12,JacPDstandardNth2gt12,kmadd(gtu22,JacPDstandardNth2gt22,kmadd(gtu23,JacPDstandardNth2gt23,kmadd(gtu13,JacPDstandardNth3gt12,kmadd(gtu23,JacPDstandardNth3gt22,kmul(gtu33,JacPDstandardNth3gt23))))))))),ToReal(-2),Ddetgt2),kmul(JacPDstandardNth2alpha,ToReal(2)))),kmul(gtu33,kmadd(alphaL,kmul(JacPDstandardNth3phi,kmul(IfThen(conformalMethod 
        == 
        1,kdiv(ToReal(1),phiL),ToReal(-2)),ToReal(-2))),kmadd(alphaL,kmadd(kmadd(gtu11,JacPDstandardNth1gt13,kmadd(gtu12,JacPDstandardNth1gt23,kmadd(gtu13,JacPDstandardNth1gt33,kmadd(gtu12,JacPDstandardNth2gt13,kmadd(gtu22,JacPDstandardNth2gt23,kmadd(gtu23,JacPDstandardNth2gt33,kmadd(gtu13,JacPDstandardNth3gt13,kmadd(gtu23,JacPDstandardNth3gt23,kmul(gtu33,JacPDstandardNth3gt33))))))))),ToReal(-2),Ddetgt3),kmul(JacPDstandardNth3alpha,ToReal(2))))))))));
    }
    else
    {
      beta1rhsL = 
        kmul(theta,kmul(kmadd(B1L,ToReal(ShiftBCoeff),kmul(knmsub(beta1L,kmul(eta,ToReal(BetaDriver)),Xt1L),ToReal(1 
        - ShiftBCoeff))),ToReal(ShiftGammaCoeff)));
      
      beta2rhsL = 
        kmul(theta,kmul(kmadd(B2L,ToReal(ShiftBCoeff),kmul(knmsub(beta2L,kmul(eta,ToReal(BetaDriver)),Xt2L),ToReal(1 
        - ShiftBCoeff))),ToReal(ShiftGammaCoeff)));
      
      beta3rhsL = 
        kmul(theta,kmul(kmadd(B3L,ToReal(ShiftBCoeff),kmul(knmsub(beta3L,kmul(eta,ToReal(BetaDriver)),Xt3L),ToReal(1 
        - ShiftBCoeff))),ToReal(ShiftGammaCoeff)));
    }
    
    CCTK_REAL_VEC B1rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(knmsub(B1L,kmul(eta,ToReal(BetaDriver)),dotXt1),ToReal(ShiftBCoeff));
    
    CCTK_REAL_VEC B2rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(knmsub(B2L,kmul(eta,ToReal(BetaDriver)),dotXt2),ToReal(ShiftBCoeff));
    
    CCTK_REAL_VEC B3rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmul(knmsub(B3L,kmul(eta,ToReal(BetaDriver)),dotXt3),ToReal(ShiftBCoeff));
    
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(alpharhs[index],alpharhsL);
    vec_store_nta_partial(Arhs[index],ArhsL);
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
  CCTK_ENDLOOP3STR(ML_BSSN_DGFE_RHS1);
}

extern "C" void ML_BSSN_DGFE_RHS1(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_DGFE_RHS1_Body");
  }
  
  if (cctk_iteration % ML_BSSN_DGFE_RHS1_calc_every != ML_BSSN_DGFE_RHS1_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "grid::coordinates",
    "ML_BSSN_DGFE::ML_curv",
    "ML_BSSN_DGFE::ML_dtlapse",
    "ML_BSSN_DGFE::ML_dtlapserhs",
    "ML_BSSN_DGFE::ML_dtshift",
    "ML_BSSN_DGFE::ML_dtshiftrhs",
    "ML_BSSN_DGFE::ML_Gamma",
    "ML_BSSN_DGFE::ML_Gammarhs",
    "ML_BSSN_DGFE::ML_lapse",
    "ML_BSSN_DGFE::ML_lapserhs",
    "ML_BSSN_DGFE::ML_log_confac",
    "ML_BSSN_DGFE::ML_log_confacrhs",
    "ML_BSSN_DGFE::ML_metric",
    "ML_BSSN_DGFE::ML_metricrhs",
    "ML_BSSN_DGFE::ML_shift",
    "ML_BSSN_DGFE::ML_shiftrhs",
    "ML_BSSN_DGFE::ML_trace_curv",
    "ML_BSSN_DGFE::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_DGFE_RHS1", 18, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_DGFE_RHS1", 1, 1, 1);
      break;
    }
    
    case 4:
    {
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_DGFE_RHS1", 2, 2, 2);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  
  if (not solver) solver = new ML_BSSN_DGFE_RHS1_method(cctkGH);
  GenericFD_LoopOverInterior(cctkGH, ML_BSSN_DGFE_RHS1_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_DGFE_RHS1_Body");
  }
}
