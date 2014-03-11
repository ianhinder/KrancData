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

extern "C" void ML_BSSN_DGFE_boundary_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_DGFE_boundary_calc_every != ML_BSSN_DGFE_boundary_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_curv","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_curv.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_dtlapse","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_dtlapse.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_dtshift","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_dtshift.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_Gamma","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_Gamma.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_lapse","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_lapse.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_log_confac","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_log_confac.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_metric","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_metric.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_shift","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_shift.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_DGFE::ML_trace_curv","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_DGFE::ML_trace_curv.");
  return;
}




/* DGFE Definitions */

#include <hrscc.hh>

#define config_sdg_order      5
#define config_riemann_solver hrscc::LaxFriedrichsRS<DGFE_ML_BSSN_DGFE_boundary, false>

/* Export definitions */
#define ML_BSSN_DGFE_boundary_sdg_grid   hrscc::GNIGrid<hrscc::GLLElement<config_sdg_order> >
#define ML_BSSN_DGFE_boundary_sdg_method hrscc::SDGMethod<DGFE_ML_BSSN_DGFE_boundary, ML_BSSN_DGFE_boundary_sdg_grid, config_riemann_solver>

/*** Numerical scheme ***/

/* Configuration */
#define config_method ML_BSSN_DGFE_boundary_sdg_method

/* Export definitions */
#define ML_BSSN_DGFE_boundary_method config_method
#define ML_BSSN_DGFE_boundary_solver hrscc::CLawSolver<DGFE_ML_BSSN_DGFE_boundary, config_method>



class DGFE_ML_BSSN_DGFE_boundary;

namespace hrscc {
  template<>
  struct traits<DGFE_ML_BSSN_DGFE_boundary> {
    // All state vector variables
    enum state_t {nvars};
    enum {nequations = nvars};
    enum {nexternal = 3*nvars};
    enum {nbitmasks = 0};
    static const bool pure = false;
  };
} // namespace



class DGFE_ML_BSSN_DGFE_boundary: public hrscc::CLaw<DGFE_ML_BSSN_DGFE_boundary> {
public:
  typedef hrscc::CLaw<DGFE_ML_BSSN_DGFE_boundary> claw;
  typedef hrscc::traits<DGFE_ML_BSSN_DGFE_boundary>::state_t state_t;
  typedef hrscc::traits<DGFE_ML_BSSN_DGFE_boundary> variables_t;
  static const int nvars = variables_t::nvars;
  
  DGFE_ML_BSSN_DGFE_boundary();
  
  inline void prim_to_all(hrscc::Observer<claw> & observer) const
  {
  }
  
  template<hrscc::policy::direction_t dir>
  inline void fluxes(hrscc::Observer<claw> & observer) const
  {
    
    
    switch (dir) {
    case hrscc::policy::x: {
      break;
    }
    case hrscc::policy::y: {
      break;
    }
    case hrscc::policy::z: {
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
    }
    
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
  template<> int CLaw<DGFE_ML_BSSN_DGFE_boundary>::conserved_idx[DGFE_ML_BSSN_DGFE_boundary::nvars] = {};
  template<> int CLaw<DGFE_ML_BSSN_DGFE_boundary>::primitive_idx[DGFE_ML_BSSN_DGFE_boundary::nvars] = {};
  template<> int CLaw<DGFE_ML_BSSN_DGFE_boundary>::rhs_idx[DGFE_ML_BSSN_DGFE_boundary::nvars] = {};
  template<> int CLaw<DGFE_ML_BSSN_DGFE_boundary>::field_idx[3*DGFE_ML_BSSN_DGFE_boundary::nvars] = {};
  template<> int CLaw<DGFE_ML_BSSN_DGFE_boundary>::bitmask_idx[0] = {};
} // namespace hrscc



namespace {
  int varindex(const char* const varname)
  {
    const int vi = CCTK_VarIndex(varname);
    if (vi<0) CCTK_ERROR("Internal error");
    return vi;
  }
}

DGFE_ML_BSSN_DGFE_boundary::DGFE_ML_BSSN_DGFE_boundary()
{
  using namespace hrscc;



}



/* A solver, DGFE's equivalent of cctkGH */
static ML_BSSN_DGFE_boundary_solver *solver = NULL;



/* Call the pointwise DGFE derivative operator */
#undef PDstandardNth1
#undef PDstandardNth2
#undef PDstandardNth3
#define PDstandardNth1(u) (solver->wdiff<hrscc::policy::x>(&(u)[-index], i,j,k))
#define PDstandardNth2(u) (solver->wdiff<hrscc::policy::y>(&(u)[-index], i,j,k))
#define PDstandardNth3(u) (solver->wdiff<hrscc::policy::z>(&(u)[-index], i,j,k))



static void ML_BSSN_DGFE_boundary_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3STR(ML_BSSN_DGFE_boundary,
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
    CCTK_REAL_VEC phiL CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod == 
      1,ToReal(1),ToReal(0));
    
    CCTK_REAL_VEC gt11L CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC gt12L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC gt13L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC gt22L CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC gt23L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC gt33L CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC trKL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At11L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At12L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At13L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At22L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At23L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC At33L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Xt1L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Xt2L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC Xt3L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC alphaL CCTK_ATTRIBUTE_UNUSED = ToReal(1);
    
    CCTK_REAL_VEC AL CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC beta1L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC beta2L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC beta3L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC B1L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC B2L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    
    CCTK_REAL_VEC B3L CCTK_ATTRIBUTE_UNUSED = ToReal(0);
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(A[index],AL);
    vec_store_nta_partial(alpha[index],alphaL);
    vec_store_nta_partial(At11[index],At11L);
    vec_store_nta_partial(At12[index],At12L);
    vec_store_nta_partial(At13[index],At13L);
    vec_store_nta_partial(At22[index],At22L);
    vec_store_nta_partial(At23[index],At23L);
    vec_store_nta_partial(At33[index],At33L);
    vec_store_nta_partial(B1[index],B1L);
    vec_store_nta_partial(B2[index],B2L);
    vec_store_nta_partial(B3[index],B3L);
    vec_store_nta_partial(beta1[index],beta1L);
    vec_store_nta_partial(beta2[index],beta2L);
    vec_store_nta_partial(beta3[index],beta3L);
    vec_store_nta_partial(gt11[index],gt11L);
    vec_store_nta_partial(gt12[index],gt12L);
    vec_store_nta_partial(gt13[index],gt13L);
    vec_store_nta_partial(gt22[index],gt22L);
    vec_store_nta_partial(gt23[index],gt23L);
    vec_store_nta_partial(gt33[index],gt33L);
    vec_store_nta_partial(phi[index],phiL);
    vec_store_nta_partial(trK[index],trKL);
    vec_store_nta_partial(Xt1[index],Xt1L);
    vec_store_nta_partial(Xt2[index],Xt2L);
    vec_store_nta_partial(Xt3[index],Xt3L);
  }
  CCTK_ENDLOOP3STR(ML_BSSN_DGFE_boundary);
}
extern "C" void ML_BSSN_DGFE_boundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_DGFE_boundary_Body");
  }
  if (cctk_iteration % ML_BSSN_DGFE_boundary_calc_every != ML_BSSN_DGFE_boundary_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ML_BSSN_DGFE::ML_curv",
    "ML_BSSN_DGFE::ML_dtlapse",
    "ML_BSSN_DGFE::ML_dtshift",
    "ML_BSSN_DGFE::ML_Gamma",
    "ML_BSSN_DGFE::ML_lapse",
    "ML_BSSN_DGFE::ML_log_confac",
    "ML_BSSN_DGFE::ML_metric",
    "ML_BSSN_DGFE::ML_shift",
    "ML_BSSN_DGFE::ML_trace_curv"};
  AssertGroupStorage(cctkGH, "ML_BSSN_DGFE_boundary", 9, groups);
  
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
  
  
  if (not solver) solver = new ML_BSSN_DGFE_boundary_method(cctkGH);
  LoopOverBoundaryWithGhosts(cctkGH, ML_BSSN_DGFE_boundary_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_DGFE_boundary_Body");
  }
}

} // namespace ML_BSSN_DGFE
