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

extern "C" void ML_BSSN_NoVec_RHS1_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_NoVec_RHS1_calc_every != ML_BSSN_NoVec_RHS1_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
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

static void ML_BSSN_NoVec_RHS1_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  #pragma omp parallel
  CCTK_LOOP3(ML_BSSN_NoVec_RHS1,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL AL CCTK_ATTRIBUTE_UNUSED = A[index];
    CCTK_REAL alphaL CCTK_ATTRIBUTE_UNUSED = alpha[index];
    CCTK_REAL At11L CCTK_ATTRIBUTE_UNUSED = At11[index];
    CCTK_REAL At12L CCTK_ATTRIBUTE_UNUSED = At12[index];
    CCTK_REAL At13L CCTK_ATTRIBUTE_UNUSED = At13[index];
    CCTK_REAL At22L CCTK_ATTRIBUTE_UNUSED = At22[index];
    CCTK_REAL At23L CCTK_ATTRIBUTE_UNUSED = At23[index];
    CCTK_REAL At33L CCTK_ATTRIBUTE_UNUSED = At33[index];
    CCTK_REAL B1L CCTK_ATTRIBUTE_UNUSED = B1[index];
    CCTK_REAL B2L CCTK_ATTRIBUTE_UNUSED = B2[index];
    CCTK_REAL B3L CCTK_ATTRIBUTE_UNUSED = B3[index];
    CCTK_REAL beta1L CCTK_ATTRIBUTE_UNUSED = beta1[index];
    CCTK_REAL beta2L CCTK_ATTRIBUTE_UNUSED = beta2[index];
    CCTK_REAL beta3L CCTK_ATTRIBUTE_UNUSED = beta3[index];
    CCTK_REAL gt11L CCTK_ATTRIBUTE_UNUSED = gt11[index];
    CCTK_REAL gt12L CCTK_ATTRIBUTE_UNUSED = gt12[index];
    CCTK_REAL gt13L CCTK_ATTRIBUTE_UNUSED = gt13[index];
    CCTK_REAL gt22L CCTK_ATTRIBUTE_UNUSED = gt22[index];
    CCTK_REAL gt23L CCTK_ATTRIBUTE_UNUSED = gt23[index];
    CCTK_REAL gt33L CCTK_ATTRIBUTE_UNUSED = gt33[index];
    CCTK_REAL phiL CCTK_ATTRIBUTE_UNUSED = phi[index];
    CCTK_REAL rL CCTK_ATTRIBUTE_UNUSED = r[index];
    CCTK_REAL trKL CCTK_ATTRIBUTE_UNUSED = trK[index];
    CCTK_REAL Xt1L CCTK_ATTRIBUTE_UNUSED = Xt1[index];
    CCTK_REAL Xt2L CCTK_ATTRIBUTE_UNUSED = Xt2[index];
    CCTK_REAL Xt3L CCTK_ATTRIBUTE_UNUSED = Xt3[index];
    
    CCTK_REAL eTttL, eTtxL, eTtyL, eTtzL, eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL CCTK_ATTRIBUTE_UNUSED ;
    
    if (assume_stress_energy_state>=0 ? assume_stress_energy_state : *stress_energy_state)
    {
      eTttL = eTtt[index];
      eTtxL = eTtx[index];
      eTtyL = eTty[index];
      eTtzL = eTtz[index];
      eTxxL = eTxx[index];
      eTxyL = eTxy[index];
      eTxzL = eTxz[index];
      eTyyL = eTyy[index];
      eTyzL = eTyz[index];
      eTzzL = eTzz[index];
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
    
    CCTK_REAL dJ111L, dJ112L, dJ113L, dJ122L, dJ123L, dJ133L, dJ211L, dJ212L, dJ213L, dJ222L, dJ223L, dJ233L, dJ311L, dJ312L, dJ313L, dJ322L, dJ323L, dJ333L, J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L CCTK_ATTRIBUTE_UNUSED ;
    
    if (use_jacobian)
    {
      dJ111L = dJ111[index];
      dJ112L = dJ112[index];
      dJ113L = dJ113[index];
      dJ122L = dJ122[index];
      dJ123L = dJ123[index];
      dJ133L = dJ133[index];
      dJ211L = dJ211[index];
      dJ212L = dJ212[index];
      dJ213L = dJ213[index];
      dJ222L = dJ222[index];
      dJ223L = dJ223[index];
      dJ233L = dJ233[index];
      dJ311L = dJ311[index];
      dJ312L = dJ312[index];
      dJ313L = dJ313[index];
      dJ322L = dJ322[index];
      dJ323L = dJ323[index];
      dJ333L = dJ333[index];
      J11L = J11[index];
      J12L = J12[index];
      J13L = J13[index];
      J21L = J21[index];
      J22L = J22[index];
      J23L = J23[index];
      J31L = J31[index];
      J32L = J32[index];
      J33L = J33[index];
    }
    /* Include user supplied include files */
    /* Precompute derivatives */
    CCTK_REAL PDstandardNth1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3trK CCTK_ATTRIBUTE_UNUSED;
    
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
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = isgn(beta1L);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = isgn(beta2L);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = isgn(beta3L);
    
    CCTK_REAL JacPDstandardNth11alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2trK CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3beta1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3trK CCTK_ATTRIBUTE_UNUSED;
    
    if (use_jacobian)
    {
      JacPDstandardNth1alpha = J11L*PDstandardNth1alpha + 
        J21L*PDstandardNth2alpha + J31L*PDstandardNth3alpha;
      
      JacPDstandardNth1beta1 = J11L*PDstandardNth1beta1 + 
        J21L*PDstandardNth2beta1 + J31L*PDstandardNth3beta1;
      
      JacPDstandardNth1beta2 = J11L*PDstandardNth1beta2 + 
        J21L*PDstandardNth2beta2 + J31L*PDstandardNth3beta2;
      
      JacPDstandardNth1beta3 = J11L*PDstandardNth1beta3 + 
        J21L*PDstandardNth2beta3 + J31L*PDstandardNth3beta3;
      
      JacPDstandardNth1gt11 = J11L*PDstandardNth1gt11 + 
        J21L*PDstandardNth2gt11 + J31L*PDstandardNth3gt11;
      
      JacPDstandardNth1gt12 = J11L*PDstandardNth1gt12 + 
        J21L*PDstandardNth2gt12 + J31L*PDstandardNth3gt12;
      
      JacPDstandardNth1gt13 = J11L*PDstandardNth1gt13 + 
        J21L*PDstandardNth2gt13 + J31L*PDstandardNth3gt13;
      
      JacPDstandardNth1gt22 = J11L*PDstandardNth1gt22 + 
        J21L*PDstandardNth2gt22 + J31L*PDstandardNth3gt22;
      
      JacPDstandardNth1gt23 = J11L*PDstandardNth1gt23 + 
        J21L*PDstandardNth2gt23 + J31L*PDstandardNth3gt23;
      
      JacPDstandardNth1gt33 = J11L*PDstandardNth1gt33 + 
        J21L*PDstandardNth2gt33 + J31L*PDstandardNth3gt33;
      
      JacPDstandardNth1phi = J11L*PDstandardNth1phi + J21L*PDstandardNth2phi 
        + J31L*PDstandardNth3phi;
      
      JacPDstandardNth1trK = J11L*PDstandardNth1trK + J21L*PDstandardNth2trK 
        + J31L*PDstandardNth3trK;
      
      JacPDstandardNth2alpha = J12L*PDstandardNth1alpha + 
        J22L*PDstandardNth2alpha + J32L*PDstandardNth3alpha;
      
      JacPDstandardNth2beta1 = J12L*PDstandardNth1beta1 + 
        J22L*PDstandardNth2beta1 + J32L*PDstandardNth3beta1;
      
      JacPDstandardNth2beta2 = J12L*PDstandardNth1beta2 + 
        J22L*PDstandardNth2beta2 + J32L*PDstandardNth3beta2;
      
      JacPDstandardNth2beta3 = J12L*PDstandardNth1beta3 + 
        J22L*PDstandardNth2beta3 + J32L*PDstandardNth3beta3;
      
      JacPDstandardNth2gt11 = J12L*PDstandardNth1gt11 + 
        J22L*PDstandardNth2gt11 + J32L*PDstandardNth3gt11;
      
      JacPDstandardNth2gt12 = J12L*PDstandardNth1gt12 + 
        J22L*PDstandardNth2gt12 + J32L*PDstandardNth3gt12;
      
      JacPDstandardNth2gt13 = J12L*PDstandardNth1gt13 + 
        J22L*PDstandardNth2gt13 + J32L*PDstandardNth3gt13;
      
      JacPDstandardNth2gt22 = J12L*PDstandardNth1gt22 + 
        J22L*PDstandardNth2gt22 + J32L*PDstandardNth3gt22;
      
      JacPDstandardNth2gt23 = J12L*PDstandardNth1gt23 + 
        J22L*PDstandardNth2gt23 + J32L*PDstandardNth3gt23;
      
      JacPDstandardNth2gt33 = J12L*PDstandardNth1gt33 + 
        J22L*PDstandardNth2gt33 + J32L*PDstandardNth3gt33;
      
      JacPDstandardNth2phi = J12L*PDstandardNth1phi + J22L*PDstandardNth2phi 
        + J32L*PDstandardNth3phi;
      
      JacPDstandardNth2trK = J12L*PDstandardNth1trK + J22L*PDstandardNth2trK 
        + J32L*PDstandardNth3trK;
      
      JacPDstandardNth3alpha = J13L*PDstandardNth1alpha + 
        J23L*PDstandardNth2alpha + J33L*PDstandardNth3alpha;
      
      JacPDstandardNth3beta1 = J13L*PDstandardNth1beta1 + 
        J23L*PDstandardNth2beta1 + J33L*PDstandardNth3beta1;
      
      JacPDstandardNth3beta2 = J13L*PDstandardNth1beta2 + 
        J23L*PDstandardNth2beta2 + J33L*PDstandardNth3beta2;
      
      JacPDstandardNth3beta3 = J13L*PDstandardNth1beta3 + 
        J23L*PDstandardNth2beta3 + J33L*PDstandardNth3beta3;
      
      JacPDstandardNth3gt11 = J13L*PDstandardNth1gt11 + 
        J23L*PDstandardNth2gt11 + J33L*PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = J13L*PDstandardNth1gt12 + 
        J23L*PDstandardNth2gt12 + J33L*PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = J13L*PDstandardNth1gt13 + 
        J23L*PDstandardNth2gt13 + J33L*PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = J13L*PDstandardNth1gt22 + 
        J23L*PDstandardNth2gt22 + J33L*PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = J13L*PDstandardNth1gt23 + 
        J23L*PDstandardNth2gt23 + J33L*PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = J13L*PDstandardNth1gt33 + 
        J23L*PDstandardNth2gt33 + J33L*PDstandardNth3gt33;
      
      JacPDstandardNth3phi = J13L*PDstandardNth1phi + J23L*PDstandardNth2phi 
        + J33L*PDstandardNth3phi;
      
      JacPDstandardNth3trK = J13L*PDstandardNth1trK + J23L*PDstandardNth2trK 
        + J33L*PDstandardNth3trK;
      
      JacPDstandardNth11alpha = 2*J11L*J21L*PDstandardNth12alpha + 
        2*J11L*J31L*PDstandardNth13alpha + dJ111L*PDstandardNth1alpha + 
        2*J21L*J31L*PDstandardNth23alpha + dJ211L*PDstandardNth2alpha + 
        dJ311L*PDstandardNth3alpha + PDstandardNth11alpha*SQR(J11L) + 
        PDstandardNth22alpha*SQR(J21L) + PDstandardNth33alpha*SQR(J31L);
      
      JacPDstandardNth11beta1 = 2*J11L*J21L*PDstandardNth12beta1 + 
        2*J11L*J31L*PDstandardNth13beta1 + dJ111L*PDstandardNth1beta1 + 
        2*J21L*J31L*PDstandardNth23beta1 + dJ211L*PDstandardNth2beta1 + 
        dJ311L*PDstandardNth3beta1 + PDstandardNth11beta1*SQR(J11L) + 
        PDstandardNth22beta1*SQR(J21L) + PDstandardNth33beta1*SQR(J31L);
      
      JacPDstandardNth11beta2 = 2*J11L*J21L*PDstandardNth12beta2 + 
        2*J11L*J31L*PDstandardNth13beta2 + dJ111L*PDstandardNth1beta2 + 
        2*J21L*J31L*PDstandardNth23beta2 + dJ211L*PDstandardNth2beta2 + 
        dJ311L*PDstandardNth3beta2 + PDstandardNth11beta2*SQR(J11L) + 
        PDstandardNth22beta2*SQR(J21L) + PDstandardNth33beta2*SQR(J31L);
      
      JacPDstandardNth11beta3 = 2*J11L*J21L*PDstandardNth12beta3 + 
        2*J11L*J31L*PDstandardNth13beta3 + dJ111L*PDstandardNth1beta3 + 
        2*J21L*J31L*PDstandardNth23beta3 + dJ211L*PDstandardNth2beta3 + 
        dJ311L*PDstandardNth3beta3 + PDstandardNth11beta3*SQR(J11L) + 
        PDstandardNth22beta3*SQR(J21L) + PDstandardNth33beta3*SQR(J31L);
      
      JacPDstandardNth22alpha = 2*J12L*J22L*PDstandardNth12alpha + 
        2*J12L*J32L*PDstandardNth13alpha + dJ122L*PDstandardNth1alpha + 
        2*J22L*J32L*PDstandardNth23alpha + dJ222L*PDstandardNth2alpha + 
        dJ322L*PDstandardNth3alpha + PDstandardNth11alpha*SQR(J12L) + 
        PDstandardNth22alpha*SQR(J22L) + PDstandardNth33alpha*SQR(J32L);
      
      JacPDstandardNth22beta1 = 2*J12L*J22L*PDstandardNth12beta1 + 
        2*J12L*J32L*PDstandardNth13beta1 + dJ122L*PDstandardNth1beta1 + 
        2*J22L*J32L*PDstandardNth23beta1 + dJ222L*PDstandardNth2beta1 + 
        dJ322L*PDstandardNth3beta1 + PDstandardNth11beta1*SQR(J12L) + 
        PDstandardNth22beta1*SQR(J22L) + PDstandardNth33beta1*SQR(J32L);
      
      JacPDstandardNth22beta2 = 2*J12L*J22L*PDstandardNth12beta2 + 
        2*J12L*J32L*PDstandardNth13beta2 + dJ122L*PDstandardNth1beta2 + 
        2*J22L*J32L*PDstandardNth23beta2 + dJ222L*PDstandardNth2beta2 + 
        dJ322L*PDstandardNth3beta2 + PDstandardNth11beta2*SQR(J12L) + 
        PDstandardNth22beta2*SQR(J22L) + PDstandardNth33beta2*SQR(J32L);
      
      JacPDstandardNth22beta3 = 2*J12L*J22L*PDstandardNth12beta3 + 
        2*J12L*J32L*PDstandardNth13beta3 + dJ122L*PDstandardNth1beta3 + 
        2*J22L*J32L*PDstandardNth23beta3 + dJ222L*PDstandardNth2beta3 + 
        dJ322L*PDstandardNth3beta3 + PDstandardNth11beta3*SQR(J12L) + 
        PDstandardNth22beta3*SQR(J22L) + PDstandardNth33beta3*SQR(J32L);
      
      JacPDstandardNth33alpha = 2*J13L*J23L*PDstandardNth12alpha + 
        2*J13L*J33L*PDstandardNth13alpha + dJ133L*PDstandardNth1alpha + 
        2*J23L*J33L*PDstandardNth23alpha + dJ233L*PDstandardNth2alpha + 
        dJ333L*PDstandardNth3alpha + PDstandardNth11alpha*SQR(J13L) + 
        PDstandardNth22alpha*SQR(J23L) + PDstandardNth33alpha*SQR(J33L);
      
      JacPDstandardNth33beta1 = 2*J13L*J23L*PDstandardNth12beta1 + 
        2*J13L*J33L*PDstandardNth13beta1 + dJ133L*PDstandardNth1beta1 + 
        2*J23L*J33L*PDstandardNth23beta1 + dJ233L*PDstandardNth2beta1 + 
        dJ333L*PDstandardNth3beta1 + PDstandardNth11beta1*SQR(J13L) + 
        PDstandardNth22beta1*SQR(J23L) + PDstandardNth33beta1*SQR(J33L);
      
      JacPDstandardNth33beta2 = 2*J13L*J23L*PDstandardNth12beta2 + 
        2*J13L*J33L*PDstandardNth13beta2 + dJ133L*PDstandardNth1beta2 + 
        2*J23L*J33L*PDstandardNth23beta2 + dJ233L*PDstandardNth2beta2 + 
        dJ333L*PDstandardNth3beta2 + PDstandardNth11beta2*SQR(J13L) + 
        PDstandardNth22beta2*SQR(J23L) + PDstandardNth33beta2*SQR(J33L);
      
      JacPDstandardNth33beta3 = 2*J13L*J23L*PDstandardNth12beta3 + 
        2*J13L*J33L*PDstandardNth13beta3 + dJ133L*PDstandardNth1beta3 + 
        2*J23L*J33L*PDstandardNth23beta3 + dJ233L*PDstandardNth2beta3 + 
        dJ333L*PDstandardNth3beta3 + PDstandardNth11beta3*SQR(J13L) + 
        PDstandardNth22beta3*SQR(J23L) + PDstandardNth33beta3*SQR(J33L);
      
      JacPDstandardNth12alpha = J11L*J12L*PDstandardNth11alpha + 
        J12L*J21L*PDstandardNth12alpha + J11L*J22L*PDstandardNth12alpha + 
        J12L*J31L*PDstandardNth13alpha + J11L*J32L*PDstandardNth13alpha + 
        dJ112L*PDstandardNth1alpha + J21L*J22L*PDstandardNth22alpha + 
        J22L*J31L*PDstandardNth23alpha + J21L*J32L*PDstandardNth23alpha + 
        dJ212L*PDstandardNth2alpha + J31L*J32L*PDstandardNth33alpha + 
        dJ312L*PDstandardNth3alpha;
      
      JacPDstandardNth12beta1 = J11L*J12L*PDstandardNth11beta1 + 
        J12L*J21L*PDstandardNth12beta1 + J11L*J22L*PDstandardNth12beta1 + 
        J12L*J31L*PDstandardNth13beta1 + J11L*J32L*PDstandardNth13beta1 + 
        dJ112L*PDstandardNth1beta1 + J21L*J22L*PDstandardNth22beta1 + 
        J22L*J31L*PDstandardNth23beta1 + J21L*J32L*PDstandardNth23beta1 + 
        dJ212L*PDstandardNth2beta1 + J31L*J32L*PDstandardNth33beta1 + 
        dJ312L*PDstandardNth3beta1;
      
      JacPDstandardNth12beta2 = J11L*J12L*PDstandardNth11beta2 + 
        J12L*J21L*PDstandardNth12beta2 + J11L*J22L*PDstandardNth12beta2 + 
        J12L*J31L*PDstandardNth13beta2 + J11L*J32L*PDstandardNth13beta2 + 
        dJ112L*PDstandardNth1beta2 + J21L*J22L*PDstandardNth22beta2 + 
        J22L*J31L*PDstandardNth23beta2 + J21L*J32L*PDstandardNth23beta2 + 
        dJ212L*PDstandardNth2beta2 + J31L*J32L*PDstandardNth33beta2 + 
        dJ312L*PDstandardNth3beta2;
      
      JacPDstandardNth12beta3 = J11L*J12L*PDstandardNth11beta3 + 
        J12L*J21L*PDstandardNth12beta3 + J11L*J22L*PDstandardNth12beta3 + 
        J12L*J31L*PDstandardNth13beta3 + J11L*J32L*PDstandardNth13beta3 + 
        dJ112L*PDstandardNth1beta3 + J21L*J22L*PDstandardNth22beta3 + 
        J22L*J31L*PDstandardNth23beta3 + J21L*J32L*PDstandardNth23beta3 + 
        dJ212L*PDstandardNth2beta3 + J31L*J32L*PDstandardNth33beta3 + 
        dJ312L*PDstandardNth3beta3;
      
      JacPDstandardNth13alpha = J11L*J13L*PDstandardNth11alpha + 
        J13L*J21L*PDstandardNth12alpha + J11L*J23L*PDstandardNth12alpha + 
        J13L*J31L*PDstandardNth13alpha + J11L*J33L*PDstandardNth13alpha + 
        dJ113L*PDstandardNth1alpha + J21L*J23L*PDstandardNth22alpha + 
        J23L*J31L*PDstandardNth23alpha + J21L*J33L*PDstandardNth23alpha + 
        dJ213L*PDstandardNth2alpha + J31L*J33L*PDstandardNth33alpha + 
        dJ313L*PDstandardNth3alpha;
      
      JacPDstandardNth13beta1 = J11L*J13L*PDstandardNth11beta1 + 
        J13L*J21L*PDstandardNth12beta1 + J11L*J23L*PDstandardNth12beta1 + 
        J13L*J31L*PDstandardNth13beta1 + J11L*J33L*PDstandardNth13beta1 + 
        dJ113L*PDstandardNth1beta1 + J21L*J23L*PDstandardNth22beta1 + 
        J23L*J31L*PDstandardNth23beta1 + J21L*J33L*PDstandardNth23beta1 + 
        dJ213L*PDstandardNth2beta1 + J31L*J33L*PDstandardNth33beta1 + 
        dJ313L*PDstandardNth3beta1;
      
      JacPDstandardNth13beta2 = J11L*J13L*PDstandardNth11beta2 + 
        J13L*J21L*PDstandardNth12beta2 + J11L*J23L*PDstandardNth12beta2 + 
        J13L*J31L*PDstandardNth13beta2 + J11L*J33L*PDstandardNth13beta2 + 
        dJ113L*PDstandardNth1beta2 + J21L*J23L*PDstandardNth22beta2 + 
        J23L*J31L*PDstandardNth23beta2 + J21L*J33L*PDstandardNth23beta2 + 
        dJ213L*PDstandardNth2beta2 + J31L*J33L*PDstandardNth33beta2 + 
        dJ313L*PDstandardNth3beta2;
      
      JacPDstandardNth13beta3 = J11L*J13L*PDstandardNth11beta3 + 
        J13L*J21L*PDstandardNth12beta3 + J11L*J23L*PDstandardNth12beta3 + 
        J13L*J31L*PDstandardNth13beta3 + J11L*J33L*PDstandardNth13beta3 + 
        dJ113L*PDstandardNth1beta3 + J21L*J23L*PDstandardNth22beta3 + 
        J23L*J31L*PDstandardNth23beta3 + J21L*J33L*PDstandardNth23beta3 + 
        dJ213L*PDstandardNth2beta3 + J31L*J33L*PDstandardNth33beta3 + 
        dJ313L*PDstandardNth3beta3;
      
      JacPDstandardNth21alpha = J11L*J12L*PDstandardNth11alpha + 
        J12L*J21L*PDstandardNth12alpha + J11L*J22L*PDstandardNth12alpha + 
        J12L*J31L*PDstandardNth13alpha + J11L*J32L*PDstandardNth13alpha + 
        dJ112L*PDstandardNth1alpha + J21L*J22L*PDstandardNth22alpha + 
        J22L*J31L*PDstandardNth23alpha + J21L*J32L*PDstandardNth23alpha + 
        dJ212L*PDstandardNth2alpha + J31L*J32L*PDstandardNth33alpha + 
        dJ312L*PDstandardNth3alpha;
      
      JacPDstandardNth21beta1 = J11L*J12L*PDstandardNth11beta1 + 
        J12L*J21L*PDstandardNth12beta1 + J11L*J22L*PDstandardNth12beta1 + 
        J12L*J31L*PDstandardNth13beta1 + J11L*J32L*PDstandardNth13beta1 + 
        dJ112L*PDstandardNth1beta1 + J21L*J22L*PDstandardNth22beta1 + 
        J22L*J31L*PDstandardNth23beta1 + J21L*J32L*PDstandardNth23beta1 + 
        dJ212L*PDstandardNth2beta1 + J31L*J32L*PDstandardNth33beta1 + 
        dJ312L*PDstandardNth3beta1;
      
      JacPDstandardNth21beta2 = J11L*J12L*PDstandardNth11beta2 + 
        J12L*J21L*PDstandardNth12beta2 + J11L*J22L*PDstandardNth12beta2 + 
        J12L*J31L*PDstandardNth13beta2 + J11L*J32L*PDstandardNth13beta2 + 
        dJ112L*PDstandardNth1beta2 + J21L*J22L*PDstandardNth22beta2 + 
        J22L*J31L*PDstandardNth23beta2 + J21L*J32L*PDstandardNth23beta2 + 
        dJ212L*PDstandardNth2beta2 + J31L*J32L*PDstandardNth33beta2 + 
        dJ312L*PDstandardNth3beta2;
      
      JacPDstandardNth21beta3 = J11L*J12L*PDstandardNth11beta3 + 
        J12L*J21L*PDstandardNth12beta3 + J11L*J22L*PDstandardNth12beta3 + 
        J12L*J31L*PDstandardNth13beta3 + J11L*J32L*PDstandardNth13beta3 + 
        dJ112L*PDstandardNth1beta3 + J21L*J22L*PDstandardNth22beta3 + 
        J22L*J31L*PDstandardNth23beta3 + J21L*J32L*PDstandardNth23beta3 + 
        dJ212L*PDstandardNth2beta3 + J31L*J32L*PDstandardNth33beta3 + 
        dJ312L*PDstandardNth3beta3;
      
      JacPDstandardNth23alpha = J12L*J13L*PDstandardNth11alpha + 
        J13L*J22L*PDstandardNth12alpha + J12L*J23L*PDstandardNth12alpha + 
        J13L*J32L*PDstandardNth13alpha + J12L*J33L*PDstandardNth13alpha + 
        dJ123L*PDstandardNth1alpha + J22L*J23L*PDstandardNth22alpha + 
        J23L*J32L*PDstandardNth23alpha + J22L*J33L*PDstandardNth23alpha + 
        dJ223L*PDstandardNth2alpha + J32L*J33L*PDstandardNth33alpha + 
        dJ323L*PDstandardNth3alpha;
      
      JacPDstandardNth23beta1 = J12L*J13L*PDstandardNth11beta1 + 
        J13L*J22L*PDstandardNth12beta1 + J12L*J23L*PDstandardNth12beta1 + 
        J13L*J32L*PDstandardNth13beta1 + J12L*J33L*PDstandardNth13beta1 + 
        dJ123L*PDstandardNth1beta1 + J22L*J23L*PDstandardNth22beta1 + 
        J23L*J32L*PDstandardNth23beta1 + J22L*J33L*PDstandardNth23beta1 + 
        dJ223L*PDstandardNth2beta1 + J32L*J33L*PDstandardNth33beta1 + 
        dJ323L*PDstandardNth3beta1;
      
      JacPDstandardNth23beta2 = J12L*J13L*PDstandardNth11beta2 + 
        J13L*J22L*PDstandardNth12beta2 + J12L*J23L*PDstandardNth12beta2 + 
        J13L*J32L*PDstandardNth13beta2 + J12L*J33L*PDstandardNth13beta2 + 
        dJ123L*PDstandardNth1beta2 + J22L*J23L*PDstandardNth22beta2 + 
        J23L*J32L*PDstandardNth23beta2 + J22L*J33L*PDstandardNth23beta2 + 
        dJ223L*PDstandardNth2beta2 + J32L*J33L*PDstandardNth33beta2 + 
        dJ323L*PDstandardNth3beta2;
      
      JacPDstandardNth23beta3 = J12L*J13L*PDstandardNth11beta3 + 
        J13L*J22L*PDstandardNth12beta3 + J12L*J23L*PDstandardNth12beta3 + 
        J13L*J32L*PDstandardNth13beta3 + J12L*J33L*PDstandardNth13beta3 + 
        dJ123L*PDstandardNth1beta3 + J22L*J23L*PDstandardNth22beta3 + 
        J23L*J32L*PDstandardNth23beta3 + J22L*J33L*PDstandardNth23beta3 + 
        dJ223L*PDstandardNth2beta3 + J32L*J33L*PDstandardNth33beta3 + 
        dJ323L*PDstandardNth3beta3;
      
      JacPDstandardNth31alpha = J11L*J13L*PDstandardNth11alpha + 
        J13L*J21L*PDstandardNth12alpha + J11L*J23L*PDstandardNth12alpha + 
        J13L*J31L*PDstandardNth13alpha + J11L*J33L*PDstandardNth13alpha + 
        dJ113L*PDstandardNth1alpha + J21L*J23L*PDstandardNth22alpha + 
        J23L*J31L*PDstandardNth23alpha + J21L*J33L*PDstandardNth23alpha + 
        dJ213L*PDstandardNth2alpha + J31L*J33L*PDstandardNth33alpha + 
        dJ313L*PDstandardNth3alpha;
      
      JacPDstandardNth31beta1 = J11L*J13L*PDstandardNth11beta1 + 
        J13L*J21L*PDstandardNth12beta1 + J11L*J23L*PDstandardNth12beta1 + 
        J13L*J31L*PDstandardNth13beta1 + J11L*J33L*PDstandardNth13beta1 + 
        dJ113L*PDstandardNth1beta1 + J21L*J23L*PDstandardNth22beta1 + 
        J23L*J31L*PDstandardNth23beta1 + J21L*J33L*PDstandardNth23beta1 + 
        dJ213L*PDstandardNth2beta1 + J31L*J33L*PDstandardNth33beta1 + 
        dJ313L*PDstandardNth3beta1;
      
      JacPDstandardNth31beta2 = J11L*J13L*PDstandardNth11beta2 + 
        J13L*J21L*PDstandardNth12beta2 + J11L*J23L*PDstandardNth12beta2 + 
        J13L*J31L*PDstandardNth13beta2 + J11L*J33L*PDstandardNth13beta2 + 
        dJ113L*PDstandardNth1beta2 + J21L*J23L*PDstandardNth22beta2 + 
        J23L*J31L*PDstandardNth23beta2 + J21L*J33L*PDstandardNth23beta2 + 
        dJ213L*PDstandardNth2beta2 + J31L*J33L*PDstandardNth33beta2 + 
        dJ313L*PDstandardNth3beta2;
      
      JacPDstandardNth31beta3 = J11L*J13L*PDstandardNth11beta3 + 
        J13L*J21L*PDstandardNth12beta3 + J11L*J23L*PDstandardNth12beta3 + 
        J13L*J31L*PDstandardNth13beta3 + J11L*J33L*PDstandardNth13beta3 + 
        dJ113L*PDstandardNth1beta3 + J21L*J23L*PDstandardNth22beta3 + 
        J23L*J31L*PDstandardNth23beta3 + J21L*J33L*PDstandardNth23beta3 + 
        dJ213L*PDstandardNth2beta3 + J31L*J33L*PDstandardNth33beta3 + 
        dJ313L*PDstandardNth3beta3;
      
      JacPDstandardNth32alpha = J12L*J13L*PDstandardNth11alpha + 
        J13L*J22L*PDstandardNth12alpha + J12L*J23L*PDstandardNth12alpha + 
        J13L*J32L*PDstandardNth13alpha + J12L*J33L*PDstandardNth13alpha + 
        dJ123L*PDstandardNth1alpha + J22L*J23L*PDstandardNth22alpha + 
        J23L*J32L*PDstandardNth23alpha + J22L*J33L*PDstandardNth23alpha + 
        dJ223L*PDstandardNth2alpha + J32L*J33L*PDstandardNth33alpha + 
        dJ323L*PDstandardNth3alpha;
      
      JacPDstandardNth32beta1 = J12L*J13L*PDstandardNth11beta1 + 
        J13L*J22L*PDstandardNth12beta1 + J12L*J23L*PDstandardNth12beta1 + 
        J13L*J32L*PDstandardNth13beta1 + J12L*J33L*PDstandardNth13beta1 + 
        dJ123L*PDstandardNth1beta1 + J22L*J23L*PDstandardNth22beta1 + 
        J23L*J32L*PDstandardNth23beta1 + J22L*J33L*PDstandardNth23beta1 + 
        dJ223L*PDstandardNth2beta1 + J32L*J33L*PDstandardNth33beta1 + 
        dJ323L*PDstandardNth3beta1;
      
      JacPDstandardNth32beta2 = J12L*J13L*PDstandardNth11beta2 + 
        J13L*J22L*PDstandardNth12beta2 + J12L*J23L*PDstandardNth12beta2 + 
        J13L*J32L*PDstandardNth13beta2 + J12L*J33L*PDstandardNth13beta2 + 
        dJ123L*PDstandardNth1beta2 + J22L*J23L*PDstandardNth22beta2 + 
        J23L*J32L*PDstandardNth23beta2 + J22L*J33L*PDstandardNth23beta2 + 
        dJ223L*PDstandardNth2beta2 + J32L*J33L*PDstandardNth33beta2 + 
        dJ323L*PDstandardNth3beta2;
      
      JacPDstandardNth32beta3 = J12L*J13L*PDstandardNth11beta3 + 
        J13L*J22L*PDstandardNth12beta3 + J12L*J23L*PDstandardNth12beta3 + 
        J13L*J32L*PDstandardNth13beta3 + J12L*J33L*PDstandardNth13beta3 + 
        dJ123L*PDstandardNth1beta3 + J22L*J23L*PDstandardNth22beta3 + 
        J23L*J32L*PDstandardNth23beta3 + J22L*J33L*PDstandardNth23beta3 + 
        dJ223L*PDstandardNth2beta3 + J32L*J33L*PDstandardNth33beta3 + 
        dJ323L*PDstandardNth3beta3;
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
    
    CCTK_REAL detgt CCTK_ATTRIBUTE_UNUSED = 1;
    
    CCTK_REAL gtu11 CCTK_ATTRIBUTE_UNUSED = INV(detgt)*(gt22L*gt33L - 
      SQR(gt23L));
    
    CCTK_REAL gtu12 CCTK_ATTRIBUTE_UNUSED = (gt13L*gt23L - 
      gt12L*gt33L)*INV(detgt);
    
    CCTK_REAL gtu13 CCTK_ATTRIBUTE_UNUSED = (-(gt13L*gt22L) + 
      gt12L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu22 CCTK_ATTRIBUTE_UNUSED = INV(detgt)*(gt11L*gt33L - 
      SQR(gt13L));
    
    CCTK_REAL gtu23 CCTK_ATTRIBUTE_UNUSED = (gt12L*gt13L - 
      gt11L*gt23L)*INV(detgt);
    
    CCTK_REAL gtu33 CCTK_ATTRIBUTE_UNUSED = INV(detgt)*(gt11L*gt22L - 
      SQR(gt12L));
    
    CCTK_REAL Gtl111 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth1gt11;
    
    CCTK_REAL Gtl112 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth2gt11;
    
    CCTK_REAL Gtl113 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth3gt11;
    
    CCTK_REAL Gtl122 CCTK_ATTRIBUTE_UNUSED = 0.5*(-JacPDstandardNth1gt22 + 
      2*JacPDstandardNth2gt12);
    
    CCTK_REAL Gtl123 CCTK_ATTRIBUTE_UNUSED = 0.5*(-JacPDstandardNth1gt23 + 
      JacPDstandardNth2gt13 + JacPDstandardNth3gt12);
    
    CCTK_REAL Gtl133 CCTK_ATTRIBUTE_UNUSED = 0.5*(-JacPDstandardNth1gt33 + 
      2*JacPDstandardNth3gt13);
    
    CCTK_REAL Gtl211 CCTK_ATTRIBUTE_UNUSED = 0.5*(2*JacPDstandardNth1gt12 
      - JacPDstandardNth2gt11);
    
    CCTK_REAL Gtl212 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth1gt22;
    
    CCTK_REAL Gtl213 CCTK_ATTRIBUTE_UNUSED = 0.5*(JacPDstandardNth1gt23 - 
      JacPDstandardNth2gt13 + JacPDstandardNth3gt12);
    
    CCTK_REAL Gtl222 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth2gt22;
    
    CCTK_REAL Gtl223 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth3gt22;
    
    CCTK_REAL Gtl233 CCTK_ATTRIBUTE_UNUSED = 0.5*(-JacPDstandardNth2gt33 + 
      2*JacPDstandardNth3gt23);
    
    CCTK_REAL Gtl311 CCTK_ATTRIBUTE_UNUSED = 0.5*(2*JacPDstandardNth1gt13 
      - JacPDstandardNth3gt11);
    
    CCTK_REAL Gtl312 CCTK_ATTRIBUTE_UNUSED = 0.5*(JacPDstandardNth1gt23 + 
      JacPDstandardNth2gt13 - JacPDstandardNth3gt12);
    
    CCTK_REAL Gtl313 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth1gt33;
    
    CCTK_REAL Gtl322 CCTK_ATTRIBUTE_UNUSED = 0.5*(2*JacPDstandardNth2gt23 
      - JacPDstandardNth3gt22);
    
    CCTK_REAL Gtl323 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth2gt33;
    
    CCTK_REAL Gtl333 CCTK_ATTRIBUTE_UNUSED = 0.5*JacPDstandardNth3gt33;
    
    CCTK_REAL Gt111 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu11 + Gtl211*gtu12 + 
      Gtl311*gtu13;
    
    CCTK_REAL Gt211 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu12 + Gtl211*gtu22 + 
      Gtl311*gtu23;
    
    CCTK_REAL Gt311 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu13 + Gtl211*gtu23 + 
      Gtl311*gtu33;
    
    CCTK_REAL Gt112 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu11 + Gtl212*gtu12 + 
      Gtl312*gtu13;
    
    CCTK_REAL Gt212 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu12 + Gtl212*gtu22 + 
      Gtl312*gtu23;
    
    CCTK_REAL Gt312 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu13 + Gtl212*gtu23 + 
      Gtl312*gtu33;
    
    CCTK_REAL Gt113 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu11 + Gtl213*gtu12 + 
      Gtl313*gtu13;
    
    CCTK_REAL Gt213 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu12 + Gtl213*gtu22 + 
      Gtl313*gtu23;
    
    CCTK_REAL Gt313 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu13 + Gtl213*gtu23 + 
      Gtl313*gtu33;
    
    CCTK_REAL Gt122 CCTK_ATTRIBUTE_UNUSED = Gtl122*gtu11 + Gtl222*gtu12 + 
      Gtl322*gtu13;
    
    CCTK_REAL Gt222 CCTK_ATTRIBUTE_UNUSED = Gtl122*gtu12 + Gtl222*gtu22 + 
      Gtl322*gtu23;
    
    CCTK_REAL Gt322 CCTK_ATTRIBUTE_UNUSED = Gtl122*gtu13 + Gtl222*gtu23 + 
      Gtl322*gtu33;
    
    CCTK_REAL Gt123 CCTK_ATTRIBUTE_UNUSED = Gtl123*gtu11 + Gtl223*gtu12 + 
      Gtl323*gtu13;
    
    CCTK_REAL Gt223 CCTK_ATTRIBUTE_UNUSED = Gtl123*gtu12 + Gtl223*gtu22 + 
      Gtl323*gtu23;
    
    CCTK_REAL Gt323 CCTK_ATTRIBUTE_UNUSED = Gtl123*gtu13 + Gtl223*gtu23 + 
      Gtl323*gtu33;
    
    CCTK_REAL Gt133 CCTK_ATTRIBUTE_UNUSED = Gtl133*gtu11 + Gtl233*gtu12 + 
      Gtl333*gtu13;
    
    CCTK_REAL Gt233 CCTK_ATTRIBUTE_UNUSED = Gtl133*gtu12 + Gtl233*gtu22 + 
      Gtl333*gtu23;
    
    CCTK_REAL Gt333 CCTK_ATTRIBUTE_UNUSED = Gtl133*gtu13 + Gtl233*gtu23 + 
      Gtl333*gtu33;
    
    CCTK_REAL Xtn1 CCTK_ATTRIBUTE_UNUSED = Gt111*gtu11 + 2*Gt112*gtu12 + 
      2*Gt113*gtu13 + Gt122*gtu22 + 2*Gt123*gtu23 + Gt133*gtu33;
    
    CCTK_REAL Xtn2 CCTK_ATTRIBUTE_UNUSED = Gt211*gtu11 + 2*Gt212*gtu12 + 
      2*Gt213*gtu13 + Gt222*gtu22 + 2*Gt223*gtu23 + Gt233*gtu33;
    
    CCTK_REAL Xtn3 CCTK_ATTRIBUTE_UNUSED = Gt311*gtu11 + 2*Gt312*gtu12 + 
      2*Gt313*gtu13 + Gt322*gtu22 + 2*Gt323*gtu23 + Gt333*gtu33;
    
    CCTK_REAL e4phi CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod == 
      1,INV(SQR(phiL)),exp(4*phiL));
    
    CCTK_REAL em4phi CCTK_ATTRIBUTE_UNUSED = INV(e4phi);
    
    CCTK_REAL fac1 CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod == 
      1,-0.5*INV(phiL),1);
    
    CCTK_REAL cdphi1 CCTK_ATTRIBUTE_UNUSED = fac1*JacPDstandardNth1phi;
    
    CCTK_REAL cdphi2 CCTK_ATTRIBUTE_UNUSED = fac1*JacPDstandardNth2phi;
    
    CCTK_REAL cdphi3 CCTK_ATTRIBUTE_UNUSED = fac1*JacPDstandardNth3phi;
    
    CCTK_REAL Atm11 CCTK_ATTRIBUTE_UNUSED = At11L*gtu11 + At12L*gtu12 + 
      At13L*gtu13;
    
    CCTK_REAL Atm21 CCTK_ATTRIBUTE_UNUSED = At11L*gtu12 + At12L*gtu22 + 
      At13L*gtu23;
    
    CCTK_REAL Atm31 CCTK_ATTRIBUTE_UNUSED = At11L*gtu13 + At12L*gtu23 + 
      At13L*gtu33;
    
    CCTK_REAL Atm12 CCTK_ATTRIBUTE_UNUSED = At12L*gtu11 + At22L*gtu12 + 
      At23L*gtu13;
    
    CCTK_REAL Atm22 CCTK_ATTRIBUTE_UNUSED = At12L*gtu12 + At22L*gtu22 + 
      At23L*gtu23;
    
    CCTK_REAL Atm32 CCTK_ATTRIBUTE_UNUSED = At12L*gtu13 + At22L*gtu23 + 
      At23L*gtu33;
    
    CCTK_REAL Atm13 CCTK_ATTRIBUTE_UNUSED = At13L*gtu11 + At23L*gtu12 + 
      At33L*gtu13;
    
    CCTK_REAL Atm23 CCTK_ATTRIBUTE_UNUSED = At13L*gtu12 + At23L*gtu22 + 
      At33L*gtu23;
    
    CCTK_REAL Atm33 CCTK_ATTRIBUTE_UNUSED = At13L*gtu13 + At23L*gtu23 + 
      At33L*gtu33;
    
    CCTK_REAL Atu11 CCTK_ATTRIBUTE_UNUSED = Atm11*gtu11 + Atm12*gtu12 + 
      Atm13*gtu13;
    
    CCTK_REAL Atu12 CCTK_ATTRIBUTE_UNUSED = Atm11*gtu12 + Atm12*gtu22 + 
      Atm13*gtu23;
    
    CCTK_REAL Atu13 CCTK_ATTRIBUTE_UNUSED = Atm11*gtu13 + Atm12*gtu23 + 
      Atm13*gtu33;
    
    CCTK_REAL Atu22 CCTK_ATTRIBUTE_UNUSED = Atm21*gtu12 + Atm22*gtu22 + 
      Atm23*gtu23;
    
    CCTK_REAL Atu23 CCTK_ATTRIBUTE_UNUSED = Atm21*gtu13 + Atm22*gtu23 + 
      Atm23*gtu33;
    
    CCTK_REAL Atu33 CCTK_ATTRIBUTE_UNUSED = Atm31*gtu13 + Atm32*gtu23 + 
      Atm33*gtu33;
    
    CCTK_REAL rho CCTK_ATTRIBUTE_UNUSED = (eTttL - 2*(beta1L*eTtxL + 
      beta2L*eTtyL + beta3L*eTtzL) + beta1L*(beta1L*eTxxL + beta2L*eTxyL + 
      beta3L*eTxzL) + beta2L*(beta1L*eTxyL + beta2L*eTyyL + beta3L*eTyzL) + 
      beta3L*(beta1L*eTxzL + beta2L*eTyzL + beta3L*eTzzL))*INV(SQR(alphaL));
    
    CCTK_REAL S1 CCTK_ATTRIBUTE_UNUSED = -((eTtxL - beta1L*eTxxL - 
      beta2L*eTxyL - beta3L*eTxzL)*INV(alphaL));
    
    CCTK_REAL S2 CCTK_ATTRIBUTE_UNUSED = -((eTtyL - beta1L*eTxyL - 
      beta2L*eTyyL - beta3L*eTyzL)*INV(alphaL));
    
    CCTK_REAL S3 CCTK_ATTRIBUTE_UNUSED = -((eTtzL - beta1L*eTxzL - 
      beta2L*eTyzL - beta3L*eTzzL)*INV(alphaL));
    
    CCTK_REAL trS CCTK_ATTRIBUTE_UNUSED = em4phi*(eTxxL*gtu11 + 
      2*eTxyL*gtu12 + 2*eTxzL*gtu13 + eTyyL*gtu22 + 2*eTyzL*gtu23 + 
      eTzzL*gtu33);
    
    CCTK_REAL phirhsL CCTK_ATTRIBUTE_UNUSED = (alphaL*trKL - 
      JacPDstandardNth1beta1 - JacPDstandardNth2beta2 - 
      JacPDstandardNth3beta3)*IfThen(conformalMethod == 
      1,0.333333333333333333333333333333*phiL,-0.166666666666666666666666666667);
    
    CCTK_REAL gt11rhsL CCTK_ATTRIBUTE_UNUSED = -2*alphaL*At11L + 
      2*gt11L*JacPDstandardNth1beta1 + 2*gt12L*JacPDstandardNth1beta2 + 
      2*gt13L*JacPDstandardNth1beta3 - 
      0.666666666666666666666666666667*gt11L*(JacPDstandardNth1beta1 + 
      JacPDstandardNth2beta2 + JacPDstandardNth3beta3);
    
    CCTK_REAL gt12rhsL CCTK_ATTRIBUTE_UNUSED = -2*alphaL*At12L + 
      gt12L*JacPDstandardNth1beta1 + gt22L*JacPDstandardNth1beta2 + 
      gt23L*JacPDstandardNth1beta3 + gt11L*JacPDstandardNth2beta1 + 
      gt12L*JacPDstandardNth2beta2 + gt13L*JacPDstandardNth2beta3 - 
      0.666666666666666666666666666667*gt12L*(JacPDstandardNth1beta1 + 
      JacPDstandardNth2beta2 + JacPDstandardNth3beta3);
    
    CCTK_REAL gt13rhsL CCTK_ATTRIBUTE_UNUSED = -2*alphaL*At13L + 
      gt13L*JacPDstandardNth1beta1 + gt23L*JacPDstandardNth1beta2 + 
      gt33L*JacPDstandardNth1beta3 + gt11L*JacPDstandardNth3beta1 + 
      gt12L*JacPDstandardNth3beta2 + gt13L*JacPDstandardNth3beta3 - 
      0.666666666666666666666666666667*gt13L*(JacPDstandardNth1beta1 + 
      JacPDstandardNth2beta2 + JacPDstandardNth3beta3);
    
    CCTK_REAL gt22rhsL CCTK_ATTRIBUTE_UNUSED = -2*alphaL*At22L + 
      2*gt12L*JacPDstandardNth2beta1 + 2*gt22L*JacPDstandardNth2beta2 + 
      2*gt23L*JacPDstandardNth2beta3 - 
      0.666666666666666666666666666667*gt22L*(JacPDstandardNth1beta1 + 
      JacPDstandardNth2beta2 + JacPDstandardNth3beta3);
    
    CCTK_REAL gt23rhsL CCTK_ATTRIBUTE_UNUSED = -2*alphaL*At23L + 
      gt13L*JacPDstandardNth2beta1 + gt23L*JacPDstandardNth2beta2 + 
      gt33L*JacPDstandardNth2beta3 + gt12L*JacPDstandardNth3beta1 + 
      gt22L*JacPDstandardNth3beta2 + gt23L*JacPDstandardNth3beta3 - 
      0.666666666666666666666666666667*gt23L*(JacPDstandardNth1beta1 + 
      JacPDstandardNth2beta2 + JacPDstandardNth3beta3);
    
    CCTK_REAL gt33rhsL CCTK_ATTRIBUTE_UNUSED = -2*alphaL*At33L + 
      2*gt13L*JacPDstandardNth3beta1 + 2*gt23L*JacPDstandardNth3beta2 + 
      2*gt33L*JacPDstandardNth3beta3 - 
      0.666666666666666666666666666667*gt33L*(JacPDstandardNth1beta1 + 
      JacPDstandardNth2beta2 + JacPDstandardNth3beta3);
    
    CCTK_REAL dotXt1 CCTK_ATTRIBUTE_UNUSED = gtu11*JacPDstandardNth11beta1 
      + gtu12*JacPDstandardNth12beta1 + gtu13*JacPDstandardNth13beta1 + 
      gtu12*JacPDstandardNth21beta1 + gtu22*JacPDstandardNth22beta1 + 
      gtu23*JacPDstandardNth23beta1 + gtu13*JacPDstandardNth31beta1 + 
      gtu23*JacPDstandardNth32beta1 + gtu33*JacPDstandardNth33beta1 + 
      0.333333333333333333333333333333*(gtu11*(JacPDstandardNth11beta1 + 
      JacPDstandardNth12beta2 + JacPDstandardNth13beta3) + 
      gtu12*(JacPDstandardNth21beta1 + JacPDstandardNth22beta2 + 
      JacPDstandardNth23beta3) + gtu13*(JacPDstandardNth31beta1 + 
      JacPDstandardNth32beta2 + JacPDstandardNth33beta3)) - 
      2*(Atu11*JacPDstandardNth1alpha + Atu12*JacPDstandardNth2alpha + 
      Atu13*JacPDstandardNth3alpha) + 2*alphaL*(6*(Atu11*cdphi1 + 
      Atu12*cdphi2 + Atu13*cdphi3) + Atu11*Gt111 + 2*Atu12*Gt112 + 
      2*Atu13*Gt113 + Atu22*Gt122 + 2*Atu23*Gt123 + Atu33*Gt133 - 
      0.666666666666666666666666666667*(gtu11*JacPDstandardNth1trK + 
      gtu12*JacPDstandardNth2trK + gtu13*JacPDstandardNth3trK)) - 
      16*alphaL*Pi*(gtu11*S1 + gtu12*S2 + gtu13*S3) - 
      JacPDstandardNth1beta1*Xtn1 + 
      0.666666666666666666666666666667*(JacPDstandardNth1beta1 + 
      JacPDstandardNth2beta2 + JacPDstandardNth3beta3)*Xtn1 - 
      JacPDstandardNth2beta1*Xtn2 - JacPDstandardNth3beta1*Xtn3;
    
    CCTK_REAL dotXt2 CCTK_ATTRIBUTE_UNUSED = gtu11*JacPDstandardNth11beta2 
      + gtu12*JacPDstandardNth12beta2 + gtu13*JacPDstandardNth13beta2 + 
      gtu12*JacPDstandardNth21beta2 + gtu22*JacPDstandardNth22beta2 + 
      gtu23*JacPDstandardNth23beta2 + gtu13*JacPDstandardNth31beta2 + 
      gtu23*JacPDstandardNth32beta2 + gtu33*JacPDstandardNth33beta2 + 
      0.333333333333333333333333333333*(gtu12*(JacPDstandardNth11beta1 + 
      JacPDstandardNth12beta2 + JacPDstandardNth13beta3) + 
      gtu22*(JacPDstandardNth21beta1 + JacPDstandardNth22beta2 + 
      JacPDstandardNth23beta3) + gtu23*(JacPDstandardNth31beta1 + 
      JacPDstandardNth32beta2 + JacPDstandardNth33beta3)) - 
      2*(Atu12*JacPDstandardNth1alpha + Atu22*JacPDstandardNth2alpha + 
      Atu23*JacPDstandardNth3alpha) + 2*alphaL*(6*(Atu12*cdphi1 + 
      Atu22*cdphi2 + Atu23*cdphi3) + Atu11*Gt211 + 2*Atu12*Gt212 + 
      2*Atu13*Gt213 + Atu22*Gt222 + 2*Atu23*Gt223 + Atu33*Gt233 - 
      0.666666666666666666666666666667*(gtu12*JacPDstandardNth1trK + 
      gtu22*JacPDstandardNth2trK + gtu23*JacPDstandardNth3trK)) - 
      16*alphaL*Pi*(gtu12*S1 + gtu22*S2 + gtu23*S3) - 
      JacPDstandardNth1beta2*Xtn1 - JacPDstandardNth2beta2*Xtn2 + 
      0.666666666666666666666666666667*(JacPDstandardNth1beta1 + 
      JacPDstandardNth2beta2 + JacPDstandardNth3beta3)*Xtn2 - 
      JacPDstandardNth3beta2*Xtn3;
    
    CCTK_REAL dotXt3 CCTK_ATTRIBUTE_UNUSED = gtu11*JacPDstandardNth11beta3 
      + gtu12*JacPDstandardNth12beta3 + gtu13*JacPDstandardNth13beta3 + 
      gtu12*JacPDstandardNth21beta3 + gtu22*JacPDstandardNth22beta3 + 
      gtu23*JacPDstandardNth23beta3 + gtu13*JacPDstandardNth31beta3 + 
      gtu23*JacPDstandardNth32beta3 + gtu33*JacPDstandardNth33beta3 + 
      0.333333333333333333333333333333*(gtu13*(JacPDstandardNth11beta1 + 
      JacPDstandardNth12beta2 + JacPDstandardNth13beta3) + 
      gtu23*(JacPDstandardNth21beta1 + JacPDstandardNth22beta2 + 
      JacPDstandardNth23beta3) + gtu33*(JacPDstandardNth31beta1 + 
      JacPDstandardNth32beta2 + JacPDstandardNth33beta3)) - 
      2*(Atu13*JacPDstandardNth1alpha + Atu23*JacPDstandardNth2alpha + 
      Atu33*JacPDstandardNth3alpha) + 2*alphaL*(6*(Atu13*cdphi1 + 
      Atu23*cdphi2 + Atu33*cdphi3) + Atu11*Gt311 + 2*Atu12*Gt312 + 
      2*Atu13*Gt313 + Atu22*Gt322 + 2*Atu23*Gt323 + Atu33*Gt333 - 
      0.666666666666666666666666666667*(gtu13*JacPDstandardNth1trK + 
      gtu23*JacPDstandardNth2trK + gtu33*JacPDstandardNth3trK)) - 
      16*alphaL*Pi*(gtu13*S1 + gtu23*S2 + gtu33*S3) - 
      JacPDstandardNth1beta3*Xtn1 - JacPDstandardNth2beta3*Xtn2 - 
      JacPDstandardNth3beta3*Xtn3 + 
      0.666666666666666666666666666667*(JacPDstandardNth1beta1 + 
      JacPDstandardNth2beta2 + JacPDstandardNth3beta3)*Xtn3;
    
    CCTK_REAL Xt1rhsL CCTK_ATTRIBUTE_UNUSED = dotXt1;
    
    CCTK_REAL Xt2rhsL CCTK_ATTRIBUTE_UNUSED = dotXt2;
    
    CCTK_REAL Xt3rhsL CCTK_ATTRIBUTE_UNUSED = dotXt3;
    
    CCTK_REAL dottrK CCTK_ATTRIBUTE_UNUSED = 4*alphaL*Pi*(rho + trS) - 
      em4phi*(gtu11*(JacPDstandardNth11alpha + 
      2*cdphi1*JacPDstandardNth1alpha) + 
      gtu12*(2*cdphi2*JacPDstandardNth1alpha + JacPDstandardNth21alpha) + 
      gtu12*(JacPDstandardNth12alpha + 2*cdphi1*JacPDstandardNth2alpha) + 
      gtu22*(JacPDstandardNth22alpha + 2*cdphi2*JacPDstandardNth2alpha) + 
      gtu13*(2*cdphi3*JacPDstandardNth1alpha + JacPDstandardNth31alpha) + 
      gtu23*(2*cdphi3*JacPDstandardNth2alpha + JacPDstandardNth32alpha) + 
      gtu13*(JacPDstandardNth13alpha + 2*cdphi1*JacPDstandardNth3alpha) + 
      gtu23*(JacPDstandardNth23alpha + 2*cdphi2*JacPDstandardNth3alpha) + 
      gtu33*(JacPDstandardNth33alpha + 2*cdphi3*JacPDstandardNth3alpha) - 
      JacPDstandardNth1alpha*Xtn1 - JacPDstandardNth2alpha*Xtn2 - 
      JacPDstandardNth3alpha*Xtn3) + alphaL*(2*Atm12*Atm21 + 2*Atm13*Atm31 + 
      2*Atm23*Atm32 + 0.333333333333333333333333333333*SQR(trKL) + SQR(Atm11) 
      + SQR(Atm22) + SQR(Atm33));
    
    CCTK_REAL trKrhsL CCTK_ATTRIBUTE_UNUSED = dottrK;
    
    CCTK_REAL alpharhsL CCTK_ATTRIBUTE_UNUSED = 
      -(pow(alphaL,ToReal(harmonicN))*ToReal(harmonicF)*((trKL + (-1 + 
      alphaL)*ToReal(AlphaDriver))*(1 - ToReal(LapseACoeff)) + 
      AL*ToReal(LapseACoeff)));
    
    CCTK_REAL ArhsL CCTK_ATTRIBUTE_UNUSED = (dottrK - 
      AL*ToReal(AlphaDriver))*ToReal(LapseACoeff);
    
    CCTK_REAL eta CCTK_ATTRIBUTE_UNUSED = 
      INV(fmax(rL,ToReal(SpatialBetaDriverRadius)))*ToReal(SpatialBetaDriverRadius);
    
    CCTK_REAL theta CCTK_ATTRIBUTE_UNUSED = fmin(1,exp(1 - 
      rL*INV(ToReal(SpatialShiftGammaCoeffRadius))));
    
    CCTK_REAL Ddetgt1 CCTK_ATTRIBUTE_UNUSED = gtu11*JacPDstandardNth1gt11 
      + 2*gtu12*JacPDstandardNth1gt12 + 2*gtu13*JacPDstandardNth1gt13 + 
      gtu22*JacPDstandardNth1gt22 + 2*gtu23*JacPDstandardNth1gt23 + 
      gtu33*JacPDstandardNth1gt33;
    
    CCTK_REAL Ddetgt2 CCTK_ATTRIBUTE_UNUSED = gtu11*JacPDstandardNth2gt11 
      + 2*gtu12*JacPDstandardNth2gt12 + 2*gtu13*JacPDstandardNth2gt13 + 
      gtu22*JacPDstandardNth2gt22 + 2*gtu23*JacPDstandardNth2gt23 + 
      gtu33*JacPDstandardNth2gt33;
    
    CCTK_REAL Ddetgt3 CCTK_ATTRIBUTE_UNUSED = gtu11*JacPDstandardNth3gt11 
      + 2*gtu12*JacPDstandardNth3gt12 + 2*gtu13*JacPDstandardNth3gt13 + 
      gtu22*JacPDstandardNth3gt22 + 2*gtu23*JacPDstandardNth3gt23 + 
      gtu33*JacPDstandardNth3gt33;
    
    CCTK_REAL beta1rhsL CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL beta2rhsL CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL beta3rhsL CCTK_ATTRIBUTE_UNUSED;
    
    if (harmonicShift)
    {
      beta1rhsL = -0.5*alphaL*em4phi*(gtu11*(2*JacPDstandardNth1alpha + 
        alphaL*(Ddetgt1 - 2*(gtu11*JacPDstandardNth1gt11 + 
        gtu12*JacPDstandardNth1gt12 + gtu13*JacPDstandardNth1gt13 + 
        gtu12*JacPDstandardNth2gt11 + gtu22*JacPDstandardNth2gt12 + 
        gtu23*JacPDstandardNth2gt13 + gtu13*JacPDstandardNth3gt11 + 
        gtu23*JacPDstandardNth3gt12 + gtu33*JacPDstandardNth3gt13)) - 
        2*alphaL*JacPDstandardNth1phi*IfThen(ToReal(conformalMethod) == 
        1,INV(phiL),-2)) + gtu12*(2*JacPDstandardNth2alpha + alphaL*(Ddetgt2 - 
        2*(gtu11*JacPDstandardNth1gt12 + gtu12*JacPDstandardNth1gt22 + 
        gtu13*JacPDstandardNth1gt23 + gtu12*JacPDstandardNth2gt12 + 
        gtu22*JacPDstandardNth2gt22 + gtu23*JacPDstandardNth2gt23 + 
        gtu13*JacPDstandardNth3gt12 + gtu23*JacPDstandardNth3gt22 + 
        gtu33*JacPDstandardNth3gt23)) - 
        2*alphaL*JacPDstandardNth2phi*IfThen(ToReal(conformalMethod) == 
        1,INV(phiL),-2)) + gtu13*(2*JacPDstandardNth3alpha + alphaL*(Ddetgt3 - 
        2*(gtu11*JacPDstandardNth1gt13 + gtu12*JacPDstandardNth1gt23 + 
        gtu13*JacPDstandardNth1gt33 + gtu12*JacPDstandardNth2gt13 + 
        gtu22*JacPDstandardNth2gt23 + gtu23*JacPDstandardNth2gt33 + 
        gtu13*JacPDstandardNth3gt13 + gtu23*JacPDstandardNth3gt23 + 
        gtu33*JacPDstandardNth3gt33)) - 
        2*alphaL*JacPDstandardNth3phi*IfThen(ToReal(conformalMethod) == 
        1,INV(phiL),-2)));
      
      beta2rhsL = -0.5*alphaL*em4phi*(gtu12*(2*JacPDstandardNth1alpha + 
        alphaL*(Ddetgt1 - 2*(gtu11*JacPDstandardNth1gt11 + 
        gtu12*JacPDstandardNth1gt12 + gtu13*JacPDstandardNth1gt13 + 
        gtu12*JacPDstandardNth2gt11 + gtu22*JacPDstandardNth2gt12 + 
        gtu23*JacPDstandardNth2gt13 + gtu13*JacPDstandardNth3gt11 + 
        gtu23*JacPDstandardNth3gt12 + gtu33*JacPDstandardNth3gt13)) - 
        2*alphaL*JacPDstandardNth1phi*IfThen(ToReal(conformalMethod) == 
        1,INV(phiL),-2)) + gtu22*(2*JacPDstandardNth2alpha + alphaL*(Ddetgt2 - 
        2*(gtu11*JacPDstandardNth1gt12 + gtu12*JacPDstandardNth1gt22 + 
        gtu13*JacPDstandardNth1gt23 + gtu12*JacPDstandardNth2gt12 + 
        gtu22*JacPDstandardNth2gt22 + gtu23*JacPDstandardNth2gt23 + 
        gtu13*JacPDstandardNth3gt12 + gtu23*JacPDstandardNth3gt22 + 
        gtu33*JacPDstandardNth3gt23)) - 
        2*alphaL*JacPDstandardNth2phi*IfThen(ToReal(conformalMethod) == 
        1,INV(phiL),-2)) + gtu23*(2*JacPDstandardNth3alpha + alphaL*(Ddetgt3 - 
        2*(gtu11*JacPDstandardNth1gt13 + gtu12*JacPDstandardNth1gt23 + 
        gtu13*JacPDstandardNth1gt33 + gtu12*JacPDstandardNth2gt13 + 
        gtu22*JacPDstandardNth2gt23 + gtu23*JacPDstandardNth2gt33 + 
        gtu13*JacPDstandardNth3gt13 + gtu23*JacPDstandardNth3gt23 + 
        gtu33*JacPDstandardNth3gt33)) - 
        2*alphaL*JacPDstandardNth3phi*IfThen(ToReal(conformalMethod) == 
        1,INV(phiL),-2)));
      
      beta3rhsL = -0.5*alphaL*em4phi*(gtu13*(2*JacPDstandardNth1alpha + 
        alphaL*(Ddetgt1 - 2*(gtu11*JacPDstandardNth1gt11 + 
        gtu12*JacPDstandardNth1gt12 + gtu13*JacPDstandardNth1gt13 + 
        gtu12*JacPDstandardNth2gt11 + gtu22*JacPDstandardNth2gt12 + 
        gtu23*JacPDstandardNth2gt13 + gtu13*JacPDstandardNth3gt11 + 
        gtu23*JacPDstandardNth3gt12 + gtu33*JacPDstandardNth3gt13)) - 
        2*alphaL*JacPDstandardNth1phi*IfThen(ToReal(conformalMethod) == 
        1,INV(phiL),-2)) + gtu23*(2*JacPDstandardNth2alpha + alphaL*(Ddetgt2 - 
        2*(gtu11*JacPDstandardNth1gt12 + gtu12*JacPDstandardNth1gt22 + 
        gtu13*JacPDstandardNth1gt23 + gtu12*JacPDstandardNth2gt12 + 
        gtu22*JacPDstandardNth2gt22 + gtu23*JacPDstandardNth2gt23 + 
        gtu13*JacPDstandardNth3gt12 + gtu23*JacPDstandardNth3gt22 + 
        gtu33*JacPDstandardNth3gt23)) - 
        2*alphaL*JacPDstandardNth2phi*IfThen(ToReal(conformalMethod) == 
        1,INV(phiL),-2)) + gtu33*(2*JacPDstandardNth3alpha + alphaL*(Ddetgt3 - 
        2*(gtu11*JacPDstandardNth1gt13 + gtu12*JacPDstandardNth1gt23 + 
        gtu13*JacPDstandardNth1gt33 + gtu12*JacPDstandardNth2gt13 + 
        gtu22*JacPDstandardNth2gt23 + gtu23*JacPDstandardNth2gt33 + 
        gtu13*JacPDstandardNth3gt13 + gtu23*JacPDstandardNth3gt23 + 
        gtu33*JacPDstandardNth3gt33)) - 
        2*alphaL*JacPDstandardNth3phi*IfThen(ToReal(conformalMethod) == 
        1,INV(phiL),-2)));
    }
    else
    {
      beta1rhsL = theta*((Xt1L - beta1L*eta*ToReal(BetaDriver))*(1 - 
        ToReal(ShiftBCoeff)) + 
        B1L*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
      
      beta2rhsL = theta*((Xt2L - beta2L*eta*ToReal(BetaDriver))*(1 - 
        ToReal(ShiftBCoeff)) + 
        B2L*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
      
      beta3rhsL = theta*((Xt3L - beta3L*eta*ToReal(BetaDriver))*(1 - 
        ToReal(ShiftBCoeff)) + 
        B3L*ToReal(ShiftBCoeff))*ToReal(ShiftGammaCoeff);
    }
    
    CCTK_REAL B1rhsL CCTK_ATTRIBUTE_UNUSED = (dotXt1 - 
      B1L*eta*ToReal(BetaDriver))*ToReal(ShiftBCoeff);
    
    CCTK_REAL B2rhsL CCTK_ATTRIBUTE_UNUSED = (dotXt2 - 
      B2L*eta*ToReal(BetaDriver))*ToReal(ShiftBCoeff);
    
    CCTK_REAL B3rhsL CCTK_ATTRIBUTE_UNUSED = (dotXt3 - 
      B3L*eta*ToReal(BetaDriver))*ToReal(ShiftBCoeff);
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    Arhs[index] = ArhsL;
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
  CCTK_ENDLOOP3(ML_BSSN_NoVec_RHS1);
}
extern "C" void ML_BSSN_NoVec_RHS1(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_NoVec_RHS1_Body");
  }
  if (cctk_iteration % ML_BSSN_NoVec_RHS1_calc_every != ML_BSSN_NoVec_RHS1_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "grid::coordinates",
    "ML_BSSN_NoVec::ML_curv",
    "ML_BSSN_NoVec::ML_dtlapse",
    "ML_BSSN_NoVec::ML_dtlapserhs",
    "ML_BSSN_NoVec::ML_dtshift",
    "ML_BSSN_NoVec::ML_dtshiftrhs",
    "ML_BSSN_NoVec::ML_Gamma",
    "ML_BSSN_NoVec::ML_Gammarhs",
    "ML_BSSN_NoVec::ML_lapse",
    "ML_BSSN_NoVec::ML_lapserhs",
    "ML_BSSN_NoVec::ML_log_confac",
    "ML_BSSN_NoVec::ML_log_confacrhs",
    "ML_BSSN_NoVec::ML_metric",
    "ML_BSSN_NoVec::ML_metricrhs",
    "ML_BSSN_NoVec::ML_shift",
    "ML_BSSN_NoVec::ML_shiftrhs",
    "ML_BSSN_NoVec::ML_trace_curv",
    "ML_BSSN_NoVec::ML_trace_curvrhs"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_NoVec_RHS1", 18, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_NoVec_RHS1", 1, 1, 1);
      break;
    }
    
    case 4:
    {
      GenericFD_EnsureStencilFits(cctkGH, "ML_BSSN_NoVec_RHS1", 2, 2, 2);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  GenericFD_LoopOverInterior(cctkGH, ML_BSSN_NoVec_RHS1_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_NoVec_RHS1_Body");
  }
}
