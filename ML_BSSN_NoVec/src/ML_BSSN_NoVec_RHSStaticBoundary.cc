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

/* Define macros used in calculations */
#define INITVALUE (42)
#define INV(x) ((CCTK_REAL)1.0 / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * SQR(x))
#define QAD(x) (SQR(SQR(x)))

extern "C" void ML_BSSN_NoVec_RHSStaticBoundary_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_NoVec_RHSStaticBoundary_calc_every != ML_BSSN_NoVec_RHSStaticBoundary_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NoVec::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_NoVec::ML_curvrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NoVec::ML_dtlapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_NoVec::ML_dtlapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NoVec::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_NoVec::ML_dtshiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NoVec::ML_Gammarhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_NoVec::ML_Gammarhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NoVec::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_NoVec::ML_lapserhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NoVec::ML_log_confacrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_NoVec::ML_log_confacrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NoVec::ML_metricrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_NoVec::ML_metricrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NoVec::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_NoVec::ML_shiftrhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NoVec::ML_trace_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_NoVec::ML_trace_curvrhs.");
  return;
}

static void ML_BSSN_NoVec_RHSStaticBoundary_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  const CCTK_REAL dx CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(0));
  const CCTK_REAL dy CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(1));
  const CCTK_REAL dz CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(2));
  const CCTK_REAL dt CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_TIME);
  const CCTK_REAL t CCTK_ATTRIBUTE_UNUSED = ToReal(cctk_time);
  const CCTK_REAL dxi CCTK_ATTRIBUTE_UNUSED = INV(dx);
  const CCTK_REAL dyi CCTK_ATTRIBUTE_UNUSED = INV(dy);
  const CCTK_REAL dzi CCTK_ATTRIBUTE_UNUSED = INV(dz);
  const CCTK_REAL khalf CCTK_ATTRIBUTE_UNUSED = 0.5;
  const CCTK_REAL kthird CCTK_ATTRIBUTE_UNUSED = 
    0.333333333333333333333333333333;
  const CCTK_REAL ktwothird CCTK_ATTRIBUTE_UNUSED = 
    0.666666666666666666666666666667;
  const CCTK_REAL kfourthird CCTK_ATTRIBUTE_UNUSED = 
    1.33333333333333333333333333333;
  const CCTK_REAL hdxi CCTK_ATTRIBUTE_UNUSED = 0.5*dxi;
  const CCTK_REAL hdyi CCTK_ATTRIBUTE_UNUSED = 0.5*dyi;
  const CCTK_REAL hdzi CCTK_ATTRIBUTE_UNUSED = 0.5*dzi;
  
  /* Initialize predefined quantities */
  const CCTK_REAL p1o12dx CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*INV(dx);
  const CCTK_REAL p1o12dy CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*INV(dy);
  const CCTK_REAL p1o12dz CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*INV(dz);
  const CCTK_REAL p1o144dxdy CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*INV(dx*dy);
  const CCTK_REAL p1o144dxdz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*INV(dx*dz);
  const CCTK_REAL p1o144dydz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*INV(dy*dz);
  const CCTK_REAL p1o24dx CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*INV(dx);
  const CCTK_REAL p1o24dy CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*INV(dy);
  const CCTK_REAL p1o24dz CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*INV(dz);
  const CCTK_REAL p1o2dx CCTK_ATTRIBUTE_UNUSED = 0.5*INV(dx);
  const CCTK_REAL p1o2dy CCTK_ATTRIBUTE_UNUSED = 0.5*INV(dy);
  const CCTK_REAL p1o2dz CCTK_ATTRIBUTE_UNUSED = 0.5*INV(dz);
  const CCTK_REAL p1o4dx CCTK_ATTRIBUTE_UNUSED = 0.25*INV(dx);
  const CCTK_REAL p1o4dxdy CCTK_ATTRIBUTE_UNUSED = 0.25*INV(dx*dy);
  const CCTK_REAL p1o4dxdz CCTK_ATTRIBUTE_UNUSED = 0.25*INV(dx*dz);
  const CCTK_REAL p1o4dy CCTK_ATTRIBUTE_UNUSED = 0.25*INV(dy);
  const CCTK_REAL p1o4dydz CCTK_ATTRIBUTE_UNUSED = 0.25*INV(dy*dz);
  const CCTK_REAL p1o4dz CCTK_ATTRIBUTE_UNUSED = 0.25*INV(dz);
  const CCTK_REAL p1o64dx CCTK_ATTRIBUTE_UNUSED = 0.015625*INV(dx);
  const CCTK_REAL p1o64dy CCTK_ATTRIBUTE_UNUSED = 0.015625*INV(dy);
  const CCTK_REAL p1o64dz CCTK_ATTRIBUTE_UNUSED = 0.015625*INV(dz);
  const CCTK_REAL p1odx CCTK_ATTRIBUTE_UNUSED = INV(dx);
  const CCTK_REAL p1odx2 CCTK_ATTRIBUTE_UNUSED = INV(SQR(dx));
  const CCTK_REAL p1ody CCTK_ATTRIBUTE_UNUSED = INV(dy);
  const CCTK_REAL p1ody2 CCTK_ATTRIBUTE_UNUSED = INV(SQR(dy));
  const CCTK_REAL p1odz CCTK_ATTRIBUTE_UNUSED = INV(dz);
  const CCTK_REAL p1odz2 CCTK_ATTRIBUTE_UNUSED = INV(SQR(dz));
  const CCTK_REAL pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*INV(SQR(dx));
  const CCTK_REAL pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*INV(SQR(dy));
  const CCTK_REAL pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*INV(SQR(dz));
  const CCTK_REAL pm1o16dx CCTK_ATTRIBUTE_UNUSED = -0.0625*INV(dx);
  const CCTK_REAL pm1o16dy CCTK_ATTRIBUTE_UNUSED = -0.0625*INV(dy);
  const CCTK_REAL pm1o16dz CCTK_ATTRIBUTE_UNUSED = -0.0625*INV(dz);
  const CCTK_REAL pm1o2dx CCTK_ATTRIBUTE_UNUSED = -0.5*INV(dx);
  const CCTK_REAL pm1o2dy CCTK_ATTRIBUTE_UNUSED = -0.5*INV(dy);
  const CCTK_REAL pm1o2dz CCTK_ATTRIBUTE_UNUSED = -0.5*INV(dz);
  const CCTK_REAL pm1o4dx CCTK_ATTRIBUTE_UNUSED = -0.25*INV(dx);
  const CCTK_REAL pm1o4dy CCTK_ATTRIBUTE_UNUSED = -0.25*INV(dy);
  const CCTK_REAL pm1o4dz CCTK_ATTRIBUTE_UNUSED = -0.25*INV(dz);
  
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
  CCTK_LOOP3(ML_BSSN_NoVec_RHSStaticBoundary,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    // ++vec_iter_counter;
    
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
    CCTK_REAL phirhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL gt11rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL gt12rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL gt13rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL gt22rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL gt23rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL gt33rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL trKrhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL At11rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL At12rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL At13rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL At22rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL At23rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL At33rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Xt1rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Xt2rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL Xt3rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL alpharhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL ArhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL beta1rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL beta2rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL beta3rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL B1rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL B2rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    CCTK_REAL B3rhsL CCTK_ATTRIBUTE_UNUSED = 0;
    
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    Arhs[index] = ArhsL;
    At11rhs[index] = At11rhsL;
    At12rhs[index] = At12rhsL;
    At13rhs[index] = At13rhsL;
    At22rhs[index] = At22rhsL;
    At23rhs[index] = At23rhsL;
    At33rhs[index] = At33rhsL;
    B1rhs[index] = B1rhsL;
    B2rhs[index] = B2rhsL;
    B3rhs[index] = B3rhsL;
    beta1rhs[index] = beta1rhsL;
    beta2rhs[index] = beta2rhsL;
    beta3rhs[index] = beta3rhsL;
    gt11rhs[index] = gt11rhsL;
    gt12rhs[index] = gt12rhsL;
    gt13rhs[index] = gt13rhsL;
    gt22rhs[index] = gt22rhsL;
    gt23rhs[index] = gt23rhsL;
    gt33rhs[index] = gt33rhsL;
    phirhs[index] = phirhsL;
    trKrhs[index] = trKrhsL;
    Xt1rhs[index] = Xt1rhsL;
    Xt2rhs[index] = Xt2rhsL;
    Xt3rhs[index] = Xt3rhsL;
  }
  CCTK_ENDLOOP3(ML_BSSN_NoVec_RHSStaticBoundary);
}

extern "C" void ML_BSSN_NoVec_RHSStaticBoundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_NoVec_RHSStaticBoundary_Body");
  }
  
  if (cctk_iteration % ML_BSSN_NoVec_RHSStaticBoundary_calc_every != ML_BSSN_NoVec_RHSStaticBoundary_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ML_BSSN_NoVec::ML_curvrhs",
    "ML_BSSN_NoVec::ML_dtlapserhs",
    "ML_BSSN_NoVec::ML_dtshiftrhs",
    "ML_BSSN_NoVec::ML_Gammarhs",
    "ML_BSSN_NoVec::ML_lapserhs",
    "ML_BSSN_NoVec::ML_log_confacrhs",
    "ML_BSSN_NoVec::ML_metricrhs",
    "ML_BSSN_NoVec::ML_shiftrhs",
    "ML_BSSN_NoVec::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_NoVec_RHSStaticBoundary", 9, groups);
  
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
  
  GenericFD_LoopOverBoundary(cctkGH, ML_BSSN_NoVec_RHSStaticBoundary_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_NoVec_RHSStaticBoundary_Body");
  }
}
