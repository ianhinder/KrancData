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

/* Define macros used in calculations */
#define INITVALUE (42)
#define INV(x) ((CCTK_REAL)1.0 / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * SQR(x))
#define QAD(x) (SQR(SQR(x)))

namespace ML_BSSN_NoVec {

extern "C" void ML_BSSN_NoVec_constraints1_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_NoVec_constraints1_calc_every != ML_BSSN_NoVec_constraints1_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_NoVec::ML_Ham","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_NoVec::ML_Ham.");
  return;
}

static void ML_BSSN_NoVec_constraints1_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(ML_BSSN_NoVec_constraints1,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL alphaL CCTK_ATTRIBUTE_UNUSED = alpha[index];
    CCTK_REAL At11L CCTK_ATTRIBUTE_UNUSED = At11[index];
    CCTK_REAL At12L CCTK_ATTRIBUTE_UNUSED = At12[index];
    CCTK_REAL At13L CCTK_ATTRIBUTE_UNUSED = At13[index];
    CCTK_REAL At22L CCTK_ATTRIBUTE_UNUSED = At22[index];
    CCTK_REAL At23L CCTK_ATTRIBUTE_UNUSED = At23[index];
    CCTK_REAL At33L CCTK_ATTRIBUTE_UNUSED = At33[index];
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
    CCTK_REAL PDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth11phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth22phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth33phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth12phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth13phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth23phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3Xt3 CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDstandardNth1gt11 = PDstandardNthfdOrder21(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder22(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder23(&gt11[index]);
        PDstandardNth11gt11 = PDstandardNthfdOrder211(&gt11[index]);
        PDstandardNth22gt11 = PDstandardNthfdOrder222(&gt11[index]);
        PDstandardNth33gt11 = PDstandardNthfdOrder233(&gt11[index]);
        PDstandardNth12gt11 = PDstandardNthfdOrder212(&gt11[index]);
        PDstandardNth13gt11 = PDstandardNthfdOrder213(&gt11[index]);
        PDstandardNth23gt11 = PDstandardNthfdOrder223(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder21(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder22(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder23(&gt12[index]);
        PDstandardNth11gt12 = PDstandardNthfdOrder211(&gt12[index]);
        PDstandardNth22gt12 = PDstandardNthfdOrder222(&gt12[index]);
        PDstandardNth33gt12 = PDstandardNthfdOrder233(&gt12[index]);
        PDstandardNth12gt12 = PDstandardNthfdOrder212(&gt12[index]);
        PDstandardNth13gt12 = PDstandardNthfdOrder213(&gt12[index]);
        PDstandardNth23gt12 = PDstandardNthfdOrder223(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder21(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder22(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder23(&gt13[index]);
        PDstandardNth11gt13 = PDstandardNthfdOrder211(&gt13[index]);
        PDstandardNth22gt13 = PDstandardNthfdOrder222(&gt13[index]);
        PDstandardNth33gt13 = PDstandardNthfdOrder233(&gt13[index]);
        PDstandardNth12gt13 = PDstandardNthfdOrder212(&gt13[index]);
        PDstandardNth13gt13 = PDstandardNthfdOrder213(&gt13[index]);
        PDstandardNth23gt13 = PDstandardNthfdOrder223(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder21(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder22(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder23(&gt22[index]);
        PDstandardNth11gt22 = PDstandardNthfdOrder211(&gt22[index]);
        PDstandardNth22gt22 = PDstandardNthfdOrder222(&gt22[index]);
        PDstandardNth33gt22 = PDstandardNthfdOrder233(&gt22[index]);
        PDstandardNth12gt22 = PDstandardNthfdOrder212(&gt22[index]);
        PDstandardNth13gt22 = PDstandardNthfdOrder213(&gt22[index]);
        PDstandardNth23gt22 = PDstandardNthfdOrder223(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder21(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder22(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder23(&gt23[index]);
        PDstandardNth11gt23 = PDstandardNthfdOrder211(&gt23[index]);
        PDstandardNth22gt23 = PDstandardNthfdOrder222(&gt23[index]);
        PDstandardNth33gt23 = PDstandardNthfdOrder233(&gt23[index]);
        PDstandardNth12gt23 = PDstandardNthfdOrder212(&gt23[index]);
        PDstandardNth13gt23 = PDstandardNthfdOrder213(&gt23[index]);
        PDstandardNth23gt23 = PDstandardNthfdOrder223(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder21(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder22(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder23(&gt33[index]);
        PDstandardNth11gt33 = PDstandardNthfdOrder211(&gt33[index]);
        PDstandardNth22gt33 = PDstandardNthfdOrder222(&gt33[index]);
        PDstandardNth33gt33 = PDstandardNthfdOrder233(&gt33[index]);
        PDstandardNth12gt33 = PDstandardNthfdOrder212(&gt33[index]);
        PDstandardNth13gt33 = PDstandardNthfdOrder213(&gt33[index]);
        PDstandardNth23gt33 = PDstandardNthfdOrder223(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder21(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder22(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder23(&phi[index]);
        PDstandardNth11phi = PDstandardNthfdOrder211(&phi[index]);
        PDstandardNth22phi = PDstandardNthfdOrder222(&phi[index]);
        PDstandardNth33phi = PDstandardNthfdOrder233(&phi[index]);
        PDstandardNth12phi = PDstandardNthfdOrder212(&phi[index]);
        PDstandardNth13phi = PDstandardNthfdOrder213(&phi[index]);
        PDstandardNth23phi = PDstandardNthfdOrder223(&phi[index]);
        PDstandardNth1Xt1 = PDstandardNthfdOrder21(&Xt1[index]);
        PDstandardNth2Xt1 = PDstandardNthfdOrder22(&Xt1[index]);
        PDstandardNth3Xt1 = PDstandardNthfdOrder23(&Xt1[index]);
        PDstandardNth1Xt2 = PDstandardNthfdOrder21(&Xt2[index]);
        PDstandardNth2Xt2 = PDstandardNthfdOrder22(&Xt2[index]);
        PDstandardNth3Xt2 = PDstandardNthfdOrder23(&Xt2[index]);
        PDstandardNth1Xt3 = PDstandardNthfdOrder21(&Xt3[index]);
        PDstandardNth2Xt3 = PDstandardNthfdOrder22(&Xt3[index]);
        PDstandardNth3Xt3 = PDstandardNthfdOrder23(&Xt3[index]);
        break;
      }
      
      case 4:
      {
        PDstandardNth1gt11 = PDstandardNthfdOrder41(&gt11[index]);
        PDstandardNth2gt11 = PDstandardNthfdOrder42(&gt11[index]);
        PDstandardNth3gt11 = PDstandardNthfdOrder43(&gt11[index]);
        PDstandardNth11gt11 = PDstandardNthfdOrder411(&gt11[index]);
        PDstandardNth22gt11 = PDstandardNthfdOrder422(&gt11[index]);
        PDstandardNth33gt11 = PDstandardNthfdOrder433(&gt11[index]);
        PDstandardNth12gt11 = PDstandardNthfdOrder412(&gt11[index]);
        PDstandardNth13gt11 = PDstandardNthfdOrder413(&gt11[index]);
        PDstandardNth23gt11 = PDstandardNthfdOrder423(&gt11[index]);
        PDstandardNth1gt12 = PDstandardNthfdOrder41(&gt12[index]);
        PDstandardNth2gt12 = PDstandardNthfdOrder42(&gt12[index]);
        PDstandardNth3gt12 = PDstandardNthfdOrder43(&gt12[index]);
        PDstandardNth11gt12 = PDstandardNthfdOrder411(&gt12[index]);
        PDstandardNth22gt12 = PDstandardNthfdOrder422(&gt12[index]);
        PDstandardNth33gt12 = PDstandardNthfdOrder433(&gt12[index]);
        PDstandardNth12gt12 = PDstandardNthfdOrder412(&gt12[index]);
        PDstandardNth13gt12 = PDstandardNthfdOrder413(&gt12[index]);
        PDstandardNth23gt12 = PDstandardNthfdOrder423(&gt12[index]);
        PDstandardNth1gt13 = PDstandardNthfdOrder41(&gt13[index]);
        PDstandardNth2gt13 = PDstandardNthfdOrder42(&gt13[index]);
        PDstandardNth3gt13 = PDstandardNthfdOrder43(&gt13[index]);
        PDstandardNth11gt13 = PDstandardNthfdOrder411(&gt13[index]);
        PDstandardNth22gt13 = PDstandardNthfdOrder422(&gt13[index]);
        PDstandardNth33gt13 = PDstandardNthfdOrder433(&gt13[index]);
        PDstandardNth12gt13 = PDstandardNthfdOrder412(&gt13[index]);
        PDstandardNth13gt13 = PDstandardNthfdOrder413(&gt13[index]);
        PDstandardNth23gt13 = PDstandardNthfdOrder423(&gt13[index]);
        PDstandardNth1gt22 = PDstandardNthfdOrder41(&gt22[index]);
        PDstandardNth2gt22 = PDstandardNthfdOrder42(&gt22[index]);
        PDstandardNth3gt22 = PDstandardNthfdOrder43(&gt22[index]);
        PDstandardNth11gt22 = PDstandardNthfdOrder411(&gt22[index]);
        PDstandardNth22gt22 = PDstandardNthfdOrder422(&gt22[index]);
        PDstandardNth33gt22 = PDstandardNthfdOrder433(&gt22[index]);
        PDstandardNth12gt22 = PDstandardNthfdOrder412(&gt22[index]);
        PDstandardNth13gt22 = PDstandardNthfdOrder413(&gt22[index]);
        PDstandardNth23gt22 = PDstandardNthfdOrder423(&gt22[index]);
        PDstandardNth1gt23 = PDstandardNthfdOrder41(&gt23[index]);
        PDstandardNth2gt23 = PDstandardNthfdOrder42(&gt23[index]);
        PDstandardNth3gt23 = PDstandardNthfdOrder43(&gt23[index]);
        PDstandardNth11gt23 = PDstandardNthfdOrder411(&gt23[index]);
        PDstandardNth22gt23 = PDstandardNthfdOrder422(&gt23[index]);
        PDstandardNth33gt23 = PDstandardNthfdOrder433(&gt23[index]);
        PDstandardNth12gt23 = PDstandardNthfdOrder412(&gt23[index]);
        PDstandardNth13gt23 = PDstandardNthfdOrder413(&gt23[index]);
        PDstandardNth23gt23 = PDstandardNthfdOrder423(&gt23[index]);
        PDstandardNth1gt33 = PDstandardNthfdOrder41(&gt33[index]);
        PDstandardNth2gt33 = PDstandardNthfdOrder42(&gt33[index]);
        PDstandardNth3gt33 = PDstandardNthfdOrder43(&gt33[index]);
        PDstandardNth11gt33 = PDstandardNthfdOrder411(&gt33[index]);
        PDstandardNth22gt33 = PDstandardNthfdOrder422(&gt33[index]);
        PDstandardNth33gt33 = PDstandardNthfdOrder433(&gt33[index]);
        PDstandardNth12gt33 = PDstandardNthfdOrder412(&gt33[index]);
        PDstandardNth13gt33 = PDstandardNthfdOrder413(&gt33[index]);
        PDstandardNth23gt33 = PDstandardNthfdOrder423(&gt33[index]);
        PDstandardNth1phi = PDstandardNthfdOrder41(&phi[index]);
        PDstandardNth2phi = PDstandardNthfdOrder42(&phi[index]);
        PDstandardNth3phi = PDstandardNthfdOrder43(&phi[index]);
        PDstandardNth11phi = PDstandardNthfdOrder411(&phi[index]);
        PDstandardNth22phi = PDstandardNthfdOrder422(&phi[index]);
        PDstandardNth33phi = PDstandardNthfdOrder433(&phi[index]);
        PDstandardNth12phi = PDstandardNthfdOrder412(&phi[index]);
        PDstandardNth13phi = PDstandardNthfdOrder413(&phi[index]);
        PDstandardNth23phi = PDstandardNthfdOrder423(&phi[index]);
        PDstandardNth1Xt1 = PDstandardNthfdOrder41(&Xt1[index]);
        PDstandardNth2Xt1 = PDstandardNthfdOrder42(&Xt1[index]);
        PDstandardNth3Xt1 = PDstandardNthfdOrder43(&Xt1[index]);
        PDstandardNth1Xt2 = PDstandardNthfdOrder41(&Xt2[index]);
        PDstandardNth2Xt2 = PDstandardNthfdOrder42(&Xt2[index]);
        PDstandardNth3Xt2 = PDstandardNthfdOrder43(&Xt2[index]);
        PDstandardNth1Xt3 = PDstandardNthfdOrder41(&Xt3[index]);
        PDstandardNth2Xt3 = PDstandardNthfdOrder42(&Xt3[index]);
        PDstandardNth3Xt3 = PDstandardNthfdOrder43(&Xt3[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL JacPDstandardNth11gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth11phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth12phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth13phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth21gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth22phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth23phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth31gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth32gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth33phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3Xt3 CCTK_ATTRIBUTE_UNUSED;
    
    if (use_jacobian)
    {
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
      
      JacPDstandardNth1Xt1 = J11L*PDstandardNth1Xt1 + J21L*PDstandardNth2Xt1 
        + J31L*PDstandardNth3Xt1;
      
      JacPDstandardNth1Xt2 = J11L*PDstandardNth1Xt2 + J21L*PDstandardNth2Xt2 
        + J31L*PDstandardNth3Xt2;
      
      JacPDstandardNth1Xt3 = J11L*PDstandardNth1Xt3 + J21L*PDstandardNth2Xt3 
        + J31L*PDstandardNth3Xt3;
      
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
      
      JacPDstandardNth2Xt1 = J12L*PDstandardNth1Xt1 + J22L*PDstandardNth2Xt1 
        + J32L*PDstandardNth3Xt1;
      
      JacPDstandardNth2Xt2 = J12L*PDstandardNth1Xt2 + J22L*PDstandardNth2Xt2 
        + J32L*PDstandardNth3Xt2;
      
      JacPDstandardNth2Xt3 = J12L*PDstandardNth1Xt3 + J22L*PDstandardNth2Xt3 
        + J32L*PDstandardNth3Xt3;
      
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
      
      JacPDstandardNth3Xt1 = J13L*PDstandardNth1Xt1 + J23L*PDstandardNth2Xt1 
        + J33L*PDstandardNth3Xt1;
      
      JacPDstandardNth3Xt2 = J13L*PDstandardNth1Xt2 + J23L*PDstandardNth2Xt2 
        + J33L*PDstandardNth3Xt2;
      
      JacPDstandardNth3Xt3 = J13L*PDstandardNth1Xt3 + J23L*PDstandardNth2Xt3 
        + J33L*PDstandardNth3Xt3;
      
      JacPDstandardNth11gt11 = 2*J11L*J21L*PDstandardNth12gt11 + 
        2*J11L*J31L*PDstandardNth13gt11 + dJ111L*PDstandardNth1gt11 + 
        2*J21L*J31L*PDstandardNth23gt11 + dJ211L*PDstandardNth2gt11 + 
        dJ311L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J11L) + 
        PDstandardNth22gt11*SQR(J21L) + PDstandardNth33gt11*SQR(J31L);
      
      JacPDstandardNth11gt12 = 2*J11L*J21L*PDstandardNth12gt12 + 
        2*J11L*J31L*PDstandardNth13gt12 + dJ111L*PDstandardNth1gt12 + 
        2*J21L*J31L*PDstandardNth23gt12 + dJ211L*PDstandardNth2gt12 + 
        dJ311L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J11L) + 
        PDstandardNth22gt12*SQR(J21L) + PDstandardNth33gt12*SQR(J31L);
      
      JacPDstandardNth11gt13 = 2*J11L*J21L*PDstandardNth12gt13 + 
        2*J11L*J31L*PDstandardNth13gt13 + dJ111L*PDstandardNth1gt13 + 
        2*J21L*J31L*PDstandardNth23gt13 + dJ211L*PDstandardNth2gt13 + 
        dJ311L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J11L) + 
        PDstandardNth22gt13*SQR(J21L) + PDstandardNth33gt13*SQR(J31L);
      
      JacPDstandardNth11gt22 = 2*J11L*J21L*PDstandardNth12gt22 + 
        2*J11L*J31L*PDstandardNth13gt22 + dJ111L*PDstandardNth1gt22 + 
        2*J21L*J31L*PDstandardNth23gt22 + dJ211L*PDstandardNth2gt22 + 
        dJ311L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J11L) + 
        PDstandardNth22gt22*SQR(J21L) + PDstandardNth33gt22*SQR(J31L);
      
      JacPDstandardNth11gt23 = 2*J11L*J21L*PDstandardNth12gt23 + 
        2*J11L*J31L*PDstandardNth13gt23 + dJ111L*PDstandardNth1gt23 + 
        2*J21L*J31L*PDstandardNth23gt23 + dJ211L*PDstandardNth2gt23 + 
        dJ311L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J11L) + 
        PDstandardNth22gt23*SQR(J21L) + PDstandardNth33gt23*SQR(J31L);
      
      JacPDstandardNth11gt33 = 2*J11L*J21L*PDstandardNth12gt33 + 
        2*J11L*J31L*PDstandardNth13gt33 + dJ111L*PDstandardNth1gt33 + 
        2*J21L*J31L*PDstandardNth23gt33 + dJ211L*PDstandardNth2gt33 + 
        dJ311L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J11L) + 
        PDstandardNth22gt33*SQR(J21L) + PDstandardNth33gt33*SQR(J31L);
      
      JacPDstandardNth11phi = 2*J11L*J21L*PDstandardNth12phi + 
        2*J11L*J31L*PDstandardNth13phi + dJ111L*PDstandardNth1phi + 
        2*J21L*J31L*PDstandardNth23phi + dJ211L*PDstandardNth2phi + 
        dJ311L*PDstandardNth3phi + PDstandardNth11phi*SQR(J11L) + 
        PDstandardNth22phi*SQR(J21L) + PDstandardNth33phi*SQR(J31L);
      
      JacPDstandardNth22gt11 = 2*J12L*J22L*PDstandardNth12gt11 + 
        2*J12L*J32L*PDstandardNth13gt11 + dJ122L*PDstandardNth1gt11 + 
        2*J22L*J32L*PDstandardNth23gt11 + dJ222L*PDstandardNth2gt11 + 
        dJ322L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J12L) + 
        PDstandardNth22gt11*SQR(J22L) + PDstandardNth33gt11*SQR(J32L);
      
      JacPDstandardNth22gt12 = 2*J12L*J22L*PDstandardNth12gt12 + 
        2*J12L*J32L*PDstandardNth13gt12 + dJ122L*PDstandardNth1gt12 + 
        2*J22L*J32L*PDstandardNth23gt12 + dJ222L*PDstandardNth2gt12 + 
        dJ322L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J12L) + 
        PDstandardNth22gt12*SQR(J22L) + PDstandardNth33gt12*SQR(J32L);
      
      JacPDstandardNth22gt13 = 2*J12L*J22L*PDstandardNth12gt13 + 
        2*J12L*J32L*PDstandardNth13gt13 + dJ122L*PDstandardNth1gt13 + 
        2*J22L*J32L*PDstandardNth23gt13 + dJ222L*PDstandardNth2gt13 + 
        dJ322L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J12L) + 
        PDstandardNth22gt13*SQR(J22L) + PDstandardNth33gt13*SQR(J32L);
      
      JacPDstandardNth22gt22 = 2*J12L*J22L*PDstandardNth12gt22 + 
        2*J12L*J32L*PDstandardNth13gt22 + dJ122L*PDstandardNth1gt22 + 
        2*J22L*J32L*PDstandardNth23gt22 + dJ222L*PDstandardNth2gt22 + 
        dJ322L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J12L) + 
        PDstandardNth22gt22*SQR(J22L) + PDstandardNth33gt22*SQR(J32L);
      
      JacPDstandardNth22gt23 = 2*J12L*J22L*PDstandardNth12gt23 + 
        2*J12L*J32L*PDstandardNth13gt23 + dJ122L*PDstandardNth1gt23 + 
        2*J22L*J32L*PDstandardNth23gt23 + dJ222L*PDstandardNth2gt23 + 
        dJ322L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J12L) + 
        PDstandardNth22gt23*SQR(J22L) + PDstandardNth33gt23*SQR(J32L);
      
      JacPDstandardNth22gt33 = 2*J12L*J22L*PDstandardNth12gt33 + 
        2*J12L*J32L*PDstandardNth13gt33 + dJ122L*PDstandardNth1gt33 + 
        2*J22L*J32L*PDstandardNth23gt33 + dJ222L*PDstandardNth2gt33 + 
        dJ322L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J12L) + 
        PDstandardNth22gt33*SQR(J22L) + PDstandardNth33gt33*SQR(J32L);
      
      JacPDstandardNth22phi = 2*J12L*J22L*PDstandardNth12phi + 
        2*J12L*J32L*PDstandardNth13phi + dJ122L*PDstandardNth1phi + 
        2*J22L*J32L*PDstandardNth23phi + dJ222L*PDstandardNth2phi + 
        dJ322L*PDstandardNth3phi + PDstandardNth11phi*SQR(J12L) + 
        PDstandardNth22phi*SQR(J22L) + PDstandardNth33phi*SQR(J32L);
      
      JacPDstandardNth33gt11 = 2*J13L*J23L*PDstandardNth12gt11 + 
        2*J13L*J33L*PDstandardNth13gt11 + dJ133L*PDstandardNth1gt11 + 
        2*J23L*J33L*PDstandardNth23gt11 + dJ233L*PDstandardNth2gt11 + 
        dJ333L*PDstandardNth3gt11 + PDstandardNth11gt11*SQR(J13L) + 
        PDstandardNth22gt11*SQR(J23L) + PDstandardNth33gt11*SQR(J33L);
      
      JacPDstandardNth33gt12 = 2*J13L*J23L*PDstandardNth12gt12 + 
        2*J13L*J33L*PDstandardNth13gt12 + dJ133L*PDstandardNth1gt12 + 
        2*J23L*J33L*PDstandardNth23gt12 + dJ233L*PDstandardNth2gt12 + 
        dJ333L*PDstandardNth3gt12 + PDstandardNth11gt12*SQR(J13L) + 
        PDstandardNth22gt12*SQR(J23L) + PDstandardNth33gt12*SQR(J33L);
      
      JacPDstandardNth33gt13 = 2*J13L*J23L*PDstandardNth12gt13 + 
        2*J13L*J33L*PDstandardNth13gt13 + dJ133L*PDstandardNth1gt13 + 
        2*J23L*J33L*PDstandardNth23gt13 + dJ233L*PDstandardNth2gt13 + 
        dJ333L*PDstandardNth3gt13 + PDstandardNth11gt13*SQR(J13L) + 
        PDstandardNth22gt13*SQR(J23L) + PDstandardNth33gt13*SQR(J33L);
      
      JacPDstandardNth33gt22 = 2*J13L*J23L*PDstandardNth12gt22 + 
        2*J13L*J33L*PDstandardNth13gt22 + dJ133L*PDstandardNth1gt22 + 
        2*J23L*J33L*PDstandardNth23gt22 + dJ233L*PDstandardNth2gt22 + 
        dJ333L*PDstandardNth3gt22 + PDstandardNth11gt22*SQR(J13L) + 
        PDstandardNth22gt22*SQR(J23L) + PDstandardNth33gt22*SQR(J33L);
      
      JacPDstandardNth33gt23 = 2*J13L*J23L*PDstandardNth12gt23 + 
        2*J13L*J33L*PDstandardNth13gt23 + dJ133L*PDstandardNth1gt23 + 
        2*J23L*J33L*PDstandardNth23gt23 + dJ233L*PDstandardNth2gt23 + 
        dJ333L*PDstandardNth3gt23 + PDstandardNth11gt23*SQR(J13L) + 
        PDstandardNth22gt23*SQR(J23L) + PDstandardNth33gt23*SQR(J33L);
      
      JacPDstandardNth33gt33 = 2*J13L*J23L*PDstandardNth12gt33 + 
        2*J13L*J33L*PDstandardNth13gt33 + dJ133L*PDstandardNth1gt33 + 
        2*J23L*J33L*PDstandardNth23gt33 + dJ233L*PDstandardNth2gt33 + 
        dJ333L*PDstandardNth3gt33 + PDstandardNth11gt33*SQR(J13L) + 
        PDstandardNth22gt33*SQR(J23L) + PDstandardNth33gt33*SQR(J33L);
      
      JacPDstandardNth33phi = 2*J13L*J23L*PDstandardNth12phi + 
        2*J13L*J33L*PDstandardNth13phi + dJ133L*PDstandardNth1phi + 
        2*J23L*J33L*PDstandardNth23phi + dJ233L*PDstandardNth2phi + 
        dJ333L*PDstandardNth3phi + PDstandardNth11phi*SQR(J13L) + 
        PDstandardNth22phi*SQR(J23L) + PDstandardNth33phi*SQR(J33L);
      
      JacPDstandardNth12gt11 = J11L*J12L*PDstandardNth11gt11 + 
        J12L*J21L*PDstandardNth12gt11 + J11L*J22L*PDstandardNth12gt11 + 
        J12L*J31L*PDstandardNth13gt11 + J11L*J32L*PDstandardNth13gt11 + 
        dJ112L*PDstandardNth1gt11 + J21L*J22L*PDstandardNth22gt11 + 
        J22L*J31L*PDstandardNth23gt11 + J21L*J32L*PDstandardNth23gt11 + 
        dJ212L*PDstandardNth2gt11 + J31L*J32L*PDstandardNth33gt11 + 
        dJ312L*PDstandardNth3gt11;
      
      JacPDstandardNth12gt12 = J11L*J12L*PDstandardNth11gt12 + 
        J12L*J21L*PDstandardNth12gt12 + J11L*J22L*PDstandardNth12gt12 + 
        J12L*J31L*PDstandardNth13gt12 + J11L*J32L*PDstandardNth13gt12 + 
        dJ112L*PDstandardNth1gt12 + J21L*J22L*PDstandardNth22gt12 + 
        J22L*J31L*PDstandardNth23gt12 + J21L*J32L*PDstandardNth23gt12 + 
        dJ212L*PDstandardNth2gt12 + J31L*J32L*PDstandardNth33gt12 + 
        dJ312L*PDstandardNth3gt12;
      
      JacPDstandardNth12gt13 = J11L*J12L*PDstandardNth11gt13 + 
        J12L*J21L*PDstandardNth12gt13 + J11L*J22L*PDstandardNth12gt13 + 
        J12L*J31L*PDstandardNth13gt13 + J11L*J32L*PDstandardNth13gt13 + 
        dJ112L*PDstandardNth1gt13 + J21L*J22L*PDstandardNth22gt13 + 
        J22L*J31L*PDstandardNth23gt13 + J21L*J32L*PDstandardNth23gt13 + 
        dJ212L*PDstandardNth2gt13 + J31L*J32L*PDstandardNth33gt13 + 
        dJ312L*PDstandardNth3gt13;
      
      JacPDstandardNth12gt22 = J11L*J12L*PDstandardNth11gt22 + 
        J12L*J21L*PDstandardNth12gt22 + J11L*J22L*PDstandardNth12gt22 + 
        J12L*J31L*PDstandardNth13gt22 + J11L*J32L*PDstandardNth13gt22 + 
        dJ112L*PDstandardNth1gt22 + J21L*J22L*PDstandardNth22gt22 + 
        J22L*J31L*PDstandardNth23gt22 + J21L*J32L*PDstandardNth23gt22 + 
        dJ212L*PDstandardNth2gt22 + J31L*J32L*PDstandardNth33gt22 + 
        dJ312L*PDstandardNth3gt22;
      
      JacPDstandardNth12gt23 = J11L*J12L*PDstandardNth11gt23 + 
        J12L*J21L*PDstandardNth12gt23 + J11L*J22L*PDstandardNth12gt23 + 
        J12L*J31L*PDstandardNth13gt23 + J11L*J32L*PDstandardNth13gt23 + 
        dJ112L*PDstandardNth1gt23 + J21L*J22L*PDstandardNth22gt23 + 
        J22L*J31L*PDstandardNth23gt23 + J21L*J32L*PDstandardNth23gt23 + 
        dJ212L*PDstandardNth2gt23 + J31L*J32L*PDstandardNth33gt23 + 
        dJ312L*PDstandardNth3gt23;
      
      JacPDstandardNth12gt33 = J11L*J12L*PDstandardNth11gt33 + 
        J12L*J21L*PDstandardNth12gt33 + J11L*J22L*PDstandardNth12gt33 + 
        J12L*J31L*PDstandardNth13gt33 + J11L*J32L*PDstandardNth13gt33 + 
        dJ112L*PDstandardNth1gt33 + J21L*J22L*PDstandardNth22gt33 + 
        J22L*J31L*PDstandardNth23gt33 + J21L*J32L*PDstandardNth23gt33 + 
        dJ212L*PDstandardNth2gt33 + J31L*J32L*PDstandardNth33gt33 + 
        dJ312L*PDstandardNth3gt33;
      
      JacPDstandardNth12phi = J11L*J12L*PDstandardNth11phi + 
        J12L*J21L*PDstandardNth12phi + J11L*J22L*PDstandardNth12phi + 
        J12L*J31L*PDstandardNth13phi + J11L*J32L*PDstandardNth13phi + 
        dJ112L*PDstandardNth1phi + J21L*J22L*PDstandardNth22phi + 
        J22L*J31L*PDstandardNth23phi + J21L*J32L*PDstandardNth23phi + 
        dJ212L*PDstandardNth2phi + J31L*J32L*PDstandardNth33phi + 
        dJ312L*PDstandardNth3phi;
      
      JacPDstandardNth13gt11 = J11L*J13L*PDstandardNth11gt11 + 
        J13L*J21L*PDstandardNth12gt11 + J11L*J23L*PDstandardNth12gt11 + 
        J13L*J31L*PDstandardNth13gt11 + J11L*J33L*PDstandardNth13gt11 + 
        dJ113L*PDstandardNth1gt11 + J21L*J23L*PDstandardNth22gt11 + 
        J23L*J31L*PDstandardNth23gt11 + J21L*J33L*PDstandardNth23gt11 + 
        dJ213L*PDstandardNth2gt11 + J31L*J33L*PDstandardNth33gt11 + 
        dJ313L*PDstandardNth3gt11;
      
      JacPDstandardNth13gt12 = J11L*J13L*PDstandardNth11gt12 + 
        J13L*J21L*PDstandardNth12gt12 + J11L*J23L*PDstandardNth12gt12 + 
        J13L*J31L*PDstandardNth13gt12 + J11L*J33L*PDstandardNth13gt12 + 
        dJ113L*PDstandardNth1gt12 + J21L*J23L*PDstandardNth22gt12 + 
        J23L*J31L*PDstandardNth23gt12 + J21L*J33L*PDstandardNth23gt12 + 
        dJ213L*PDstandardNth2gt12 + J31L*J33L*PDstandardNth33gt12 + 
        dJ313L*PDstandardNth3gt12;
      
      JacPDstandardNth13gt13 = J11L*J13L*PDstandardNth11gt13 + 
        J13L*J21L*PDstandardNth12gt13 + J11L*J23L*PDstandardNth12gt13 + 
        J13L*J31L*PDstandardNth13gt13 + J11L*J33L*PDstandardNth13gt13 + 
        dJ113L*PDstandardNth1gt13 + J21L*J23L*PDstandardNth22gt13 + 
        J23L*J31L*PDstandardNth23gt13 + J21L*J33L*PDstandardNth23gt13 + 
        dJ213L*PDstandardNth2gt13 + J31L*J33L*PDstandardNth33gt13 + 
        dJ313L*PDstandardNth3gt13;
      
      JacPDstandardNth13gt22 = J11L*J13L*PDstandardNth11gt22 + 
        J13L*J21L*PDstandardNth12gt22 + J11L*J23L*PDstandardNth12gt22 + 
        J13L*J31L*PDstandardNth13gt22 + J11L*J33L*PDstandardNth13gt22 + 
        dJ113L*PDstandardNth1gt22 + J21L*J23L*PDstandardNth22gt22 + 
        J23L*J31L*PDstandardNth23gt22 + J21L*J33L*PDstandardNth23gt22 + 
        dJ213L*PDstandardNth2gt22 + J31L*J33L*PDstandardNth33gt22 + 
        dJ313L*PDstandardNth3gt22;
      
      JacPDstandardNth13gt23 = J11L*J13L*PDstandardNth11gt23 + 
        J13L*J21L*PDstandardNth12gt23 + J11L*J23L*PDstandardNth12gt23 + 
        J13L*J31L*PDstandardNth13gt23 + J11L*J33L*PDstandardNth13gt23 + 
        dJ113L*PDstandardNth1gt23 + J21L*J23L*PDstandardNth22gt23 + 
        J23L*J31L*PDstandardNth23gt23 + J21L*J33L*PDstandardNth23gt23 + 
        dJ213L*PDstandardNth2gt23 + J31L*J33L*PDstandardNth33gt23 + 
        dJ313L*PDstandardNth3gt23;
      
      JacPDstandardNth13gt33 = J11L*J13L*PDstandardNth11gt33 + 
        J13L*J21L*PDstandardNth12gt33 + J11L*J23L*PDstandardNth12gt33 + 
        J13L*J31L*PDstandardNth13gt33 + J11L*J33L*PDstandardNth13gt33 + 
        dJ113L*PDstandardNth1gt33 + J21L*J23L*PDstandardNth22gt33 + 
        J23L*J31L*PDstandardNth23gt33 + J21L*J33L*PDstandardNth23gt33 + 
        dJ213L*PDstandardNth2gt33 + J31L*J33L*PDstandardNth33gt33 + 
        dJ313L*PDstandardNth3gt33;
      
      JacPDstandardNth13phi = J11L*J13L*PDstandardNth11phi + 
        J13L*J21L*PDstandardNth12phi + J11L*J23L*PDstandardNth12phi + 
        J13L*J31L*PDstandardNth13phi + J11L*J33L*PDstandardNth13phi + 
        dJ113L*PDstandardNth1phi + J21L*J23L*PDstandardNth22phi + 
        J23L*J31L*PDstandardNth23phi + J21L*J33L*PDstandardNth23phi + 
        dJ213L*PDstandardNth2phi + J31L*J33L*PDstandardNth33phi + 
        dJ313L*PDstandardNth3phi;
      
      JacPDstandardNth21gt11 = J11L*J12L*PDstandardNth11gt11 + 
        J12L*J21L*PDstandardNth12gt11 + J11L*J22L*PDstandardNth12gt11 + 
        J12L*J31L*PDstandardNth13gt11 + J11L*J32L*PDstandardNth13gt11 + 
        dJ112L*PDstandardNth1gt11 + J21L*J22L*PDstandardNth22gt11 + 
        J22L*J31L*PDstandardNth23gt11 + J21L*J32L*PDstandardNth23gt11 + 
        dJ212L*PDstandardNth2gt11 + J31L*J32L*PDstandardNth33gt11 + 
        dJ312L*PDstandardNth3gt11;
      
      JacPDstandardNth21gt12 = J11L*J12L*PDstandardNth11gt12 + 
        J12L*J21L*PDstandardNth12gt12 + J11L*J22L*PDstandardNth12gt12 + 
        J12L*J31L*PDstandardNth13gt12 + J11L*J32L*PDstandardNth13gt12 + 
        dJ112L*PDstandardNth1gt12 + J21L*J22L*PDstandardNth22gt12 + 
        J22L*J31L*PDstandardNth23gt12 + J21L*J32L*PDstandardNth23gt12 + 
        dJ212L*PDstandardNth2gt12 + J31L*J32L*PDstandardNth33gt12 + 
        dJ312L*PDstandardNth3gt12;
      
      JacPDstandardNth21gt13 = J11L*J12L*PDstandardNth11gt13 + 
        J12L*J21L*PDstandardNth12gt13 + J11L*J22L*PDstandardNth12gt13 + 
        J12L*J31L*PDstandardNth13gt13 + J11L*J32L*PDstandardNth13gt13 + 
        dJ112L*PDstandardNth1gt13 + J21L*J22L*PDstandardNth22gt13 + 
        J22L*J31L*PDstandardNth23gt13 + J21L*J32L*PDstandardNth23gt13 + 
        dJ212L*PDstandardNth2gt13 + J31L*J32L*PDstandardNth33gt13 + 
        dJ312L*PDstandardNth3gt13;
      
      JacPDstandardNth21gt22 = J11L*J12L*PDstandardNth11gt22 + 
        J12L*J21L*PDstandardNth12gt22 + J11L*J22L*PDstandardNth12gt22 + 
        J12L*J31L*PDstandardNth13gt22 + J11L*J32L*PDstandardNth13gt22 + 
        dJ112L*PDstandardNth1gt22 + J21L*J22L*PDstandardNth22gt22 + 
        J22L*J31L*PDstandardNth23gt22 + J21L*J32L*PDstandardNth23gt22 + 
        dJ212L*PDstandardNth2gt22 + J31L*J32L*PDstandardNth33gt22 + 
        dJ312L*PDstandardNth3gt22;
      
      JacPDstandardNth21gt23 = J11L*J12L*PDstandardNth11gt23 + 
        J12L*J21L*PDstandardNth12gt23 + J11L*J22L*PDstandardNth12gt23 + 
        J12L*J31L*PDstandardNth13gt23 + J11L*J32L*PDstandardNth13gt23 + 
        dJ112L*PDstandardNth1gt23 + J21L*J22L*PDstandardNth22gt23 + 
        J22L*J31L*PDstandardNth23gt23 + J21L*J32L*PDstandardNth23gt23 + 
        dJ212L*PDstandardNth2gt23 + J31L*J32L*PDstandardNth33gt23 + 
        dJ312L*PDstandardNth3gt23;
      
      JacPDstandardNth21gt33 = J11L*J12L*PDstandardNth11gt33 + 
        J12L*J21L*PDstandardNth12gt33 + J11L*J22L*PDstandardNth12gt33 + 
        J12L*J31L*PDstandardNth13gt33 + J11L*J32L*PDstandardNth13gt33 + 
        dJ112L*PDstandardNth1gt33 + J21L*J22L*PDstandardNth22gt33 + 
        J22L*J31L*PDstandardNth23gt33 + J21L*J32L*PDstandardNth23gt33 + 
        dJ212L*PDstandardNth2gt33 + J31L*J32L*PDstandardNth33gt33 + 
        dJ312L*PDstandardNth3gt33;
      
      JacPDstandardNth23gt11 = J12L*J13L*PDstandardNth11gt11 + 
        J13L*J22L*PDstandardNth12gt11 + J12L*J23L*PDstandardNth12gt11 + 
        J13L*J32L*PDstandardNth13gt11 + J12L*J33L*PDstandardNth13gt11 + 
        dJ123L*PDstandardNth1gt11 + J22L*J23L*PDstandardNth22gt11 + 
        J23L*J32L*PDstandardNth23gt11 + J22L*J33L*PDstandardNth23gt11 + 
        dJ223L*PDstandardNth2gt11 + J32L*J33L*PDstandardNth33gt11 + 
        dJ323L*PDstandardNth3gt11;
      
      JacPDstandardNth23gt12 = J12L*J13L*PDstandardNth11gt12 + 
        J13L*J22L*PDstandardNth12gt12 + J12L*J23L*PDstandardNth12gt12 + 
        J13L*J32L*PDstandardNth13gt12 + J12L*J33L*PDstandardNth13gt12 + 
        dJ123L*PDstandardNth1gt12 + J22L*J23L*PDstandardNth22gt12 + 
        J23L*J32L*PDstandardNth23gt12 + J22L*J33L*PDstandardNth23gt12 + 
        dJ223L*PDstandardNth2gt12 + J32L*J33L*PDstandardNth33gt12 + 
        dJ323L*PDstandardNth3gt12;
      
      JacPDstandardNth23gt13 = J12L*J13L*PDstandardNth11gt13 + 
        J13L*J22L*PDstandardNth12gt13 + J12L*J23L*PDstandardNth12gt13 + 
        J13L*J32L*PDstandardNth13gt13 + J12L*J33L*PDstandardNth13gt13 + 
        dJ123L*PDstandardNth1gt13 + J22L*J23L*PDstandardNth22gt13 + 
        J23L*J32L*PDstandardNth23gt13 + J22L*J33L*PDstandardNth23gt13 + 
        dJ223L*PDstandardNth2gt13 + J32L*J33L*PDstandardNth33gt13 + 
        dJ323L*PDstandardNth3gt13;
      
      JacPDstandardNth23gt22 = J12L*J13L*PDstandardNth11gt22 + 
        J13L*J22L*PDstandardNth12gt22 + J12L*J23L*PDstandardNth12gt22 + 
        J13L*J32L*PDstandardNth13gt22 + J12L*J33L*PDstandardNth13gt22 + 
        dJ123L*PDstandardNth1gt22 + J22L*J23L*PDstandardNth22gt22 + 
        J23L*J32L*PDstandardNth23gt22 + J22L*J33L*PDstandardNth23gt22 + 
        dJ223L*PDstandardNth2gt22 + J32L*J33L*PDstandardNth33gt22 + 
        dJ323L*PDstandardNth3gt22;
      
      JacPDstandardNth23gt23 = J12L*J13L*PDstandardNth11gt23 + 
        J13L*J22L*PDstandardNth12gt23 + J12L*J23L*PDstandardNth12gt23 + 
        J13L*J32L*PDstandardNth13gt23 + J12L*J33L*PDstandardNth13gt23 + 
        dJ123L*PDstandardNth1gt23 + J22L*J23L*PDstandardNth22gt23 + 
        J23L*J32L*PDstandardNth23gt23 + J22L*J33L*PDstandardNth23gt23 + 
        dJ223L*PDstandardNth2gt23 + J32L*J33L*PDstandardNth33gt23 + 
        dJ323L*PDstandardNth3gt23;
      
      JacPDstandardNth23gt33 = J12L*J13L*PDstandardNth11gt33 + 
        J13L*J22L*PDstandardNth12gt33 + J12L*J23L*PDstandardNth12gt33 + 
        J13L*J32L*PDstandardNth13gt33 + J12L*J33L*PDstandardNth13gt33 + 
        dJ123L*PDstandardNth1gt33 + J22L*J23L*PDstandardNth22gt33 + 
        J23L*J32L*PDstandardNth23gt33 + J22L*J33L*PDstandardNth23gt33 + 
        dJ223L*PDstandardNth2gt33 + J32L*J33L*PDstandardNth33gt33 + 
        dJ323L*PDstandardNth3gt33;
      
      JacPDstandardNth23phi = J12L*J13L*PDstandardNth11phi + 
        J13L*J22L*PDstandardNth12phi + J12L*J23L*PDstandardNth12phi + 
        J13L*J32L*PDstandardNth13phi + J12L*J33L*PDstandardNth13phi + 
        dJ123L*PDstandardNth1phi + J22L*J23L*PDstandardNth22phi + 
        J23L*J32L*PDstandardNth23phi + J22L*J33L*PDstandardNth23phi + 
        dJ223L*PDstandardNth2phi + J32L*J33L*PDstandardNth33phi + 
        dJ323L*PDstandardNth3phi;
      
      JacPDstandardNth31gt11 = J11L*J13L*PDstandardNth11gt11 + 
        J13L*J21L*PDstandardNth12gt11 + J11L*J23L*PDstandardNth12gt11 + 
        J13L*J31L*PDstandardNth13gt11 + J11L*J33L*PDstandardNth13gt11 + 
        dJ113L*PDstandardNth1gt11 + J21L*J23L*PDstandardNth22gt11 + 
        J23L*J31L*PDstandardNth23gt11 + J21L*J33L*PDstandardNth23gt11 + 
        dJ213L*PDstandardNth2gt11 + J31L*J33L*PDstandardNth33gt11 + 
        dJ313L*PDstandardNth3gt11;
      
      JacPDstandardNth31gt12 = J11L*J13L*PDstandardNth11gt12 + 
        J13L*J21L*PDstandardNth12gt12 + J11L*J23L*PDstandardNth12gt12 + 
        J13L*J31L*PDstandardNth13gt12 + J11L*J33L*PDstandardNth13gt12 + 
        dJ113L*PDstandardNth1gt12 + J21L*J23L*PDstandardNth22gt12 + 
        J23L*J31L*PDstandardNth23gt12 + J21L*J33L*PDstandardNth23gt12 + 
        dJ213L*PDstandardNth2gt12 + J31L*J33L*PDstandardNth33gt12 + 
        dJ313L*PDstandardNth3gt12;
      
      JacPDstandardNth31gt13 = J11L*J13L*PDstandardNth11gt13 + 
        J13L*J21L*PDstandardNth12gt13 + J11L*J23L*PDstandardNth12gt13 + 
        J13L*J31L*PDstandardNth13gt13 + J11L*J33L*PDstandardNth13gt13 + 
        dJ113L*PDstandardNth1gt13 + J21L*J23L*PDstandardNth22gt13 + 
        J23L*J31L*PDstandardNth23gt13 + J21L*J33L*PDstandardNth23gt13 + 
        dJ213L*PDstandardNth2gt13 + J31L*J33L*PDstandardNth33gt13 + 
        dJ313L*PDstandardNth3gt13;
      
      JacPDstandardNth31gt22 = J11L*J13L*PDstandardNth11gt22 + 
        J13L*J21L*PDstandardNth12gt22 + J11L*J23L*PDstandardNth12gt22 + 
        J13L*J31L*PDstandardNth13gt22 + J11L*J33L*PDstandardNth13gt22 + 
        dJ113L*PDstandardNth1gt22 + J21L*J23L*PDstandardNth22gt22 + 
        J23L*J31L*PDstandardNth23gt22 + J21L*J33L*PDstandardNth23gt22 + 
        dJ213L*PDstandardNth2gt22 + J31L*J33L*PDstandardNth33gt22 + 
        dJ313L*PDstandardNth3gt22;
      
      JacPDstandardNth31gt23 = J11L*J13L*PDstandardNth11gt23 + 
        J13L*J21L*PDstandardNth12gt23 + J11L*J23L*PDstandardNth12gt23 + 
        J13L*J31L*PDstandardNth13gt23 + J11L*J33L*PDstandardNth13gt23 + 
        dJ113L*PDstandardNth1gt23 + J21L*J23L*PDstandardNth22gt23 + 
        J23L*J31L*PDstandardNth23gt23 + J21L*J33L*PDstandardNth23gt23 + 
        dJ213L*PDstandardNth2gt23 + J31L*J33L*PDstandardNth33gt23 + 
        dJ313L*PDstandardNth3gt23;
      
      JacPDstandardNth31gt33 = J11L*J13L*PDstandardNth11gt33 + 
        J13L*J21L*PDstandardNth12gt33 + J11L*J23L*PDstandardNth12gt33 + 
        J13L*J31L*PDstandardNth13gt33 + J11L*J33L*PDstandardNth13gt33 + 
        dJ113L*PDstandardNth1gt33 + J21L*J23L*PDstandardNth22gt33 + 
        J23L*J31L*PDstandardNth23gt33 + J21L*J33L*PDstandardNth23gt33 + 
        dJ213L*PDstandardNth2gt33 + J31L*J33L*PDstandardNth33gt33 + 
        dJ313L*PDstandardNth3gt33;
      
      JacPDstandardNth32gt11 = J12L*J13L*PDstandardNth11gt11 + 
        J13L*J22L*PDstandardNth12gt11 + J12L*J23L*PDstandardNth12gt11 + 
        J13L*J32L*PDstandardNth13gt11 + J12L*J33L*PDstandardNth13gt11 + 
        dJ123L*PDstandardNth1gt11 + J22L*J23L*PDstandardNth22gt11 + 
        J23L*J32L*PDstandardNth23gt11 + J22L*J33L*PDstandardNth23gt11 + 
        dJ223L*PDstandardNth2gt11 + J32L*J33L*PDstandardNth33gt11 + 
        dJ323L*PDstandardNth3gt11;
      
      JacPDstandardNth32gt12 = J12L*J13L*PDstandardNth11gt12 + 
        J13L*J22L*PDstandardNth12gt12 + J12L*J23L*PDstandardNth12gt12 + 
        J13L*J32L*PDstandardNth13gt12 + J12L*J33L*PDstandardNth13gt12 + 
        dJ123L*PDstandardNth1gt12 + J22L*J23L*PDstandardNth22gt12 + 
        J23L*J32L*PDstandardNth23gt12 + J22L*J33L*PDstandardNth23gt12 + 
        dJ223L*PDstandardNth2gt12 + J32L*J33L*PDstandardNth33gt12 + 
        dJ323L*PDstandardNth3gt12;
      
      JacPDstandardNth32gt13 = J12L*J13L*PDstandardNth11gt13 + 
        J13L*J22L*PDstandardNth12gt13 + J12L*J23L*PDstandardNth12gt13 + 
        J13L*J32L*PDstandardNth13gt13 + J12L*J33L*PDstandardNth13gt13 + 
        dJ123L*PDstandardNth1gt13 + J22L*J23L*PDstandardNth22gt13 + 
        J23L*J32L*PDstandardNth23gt13 + J22L*J33L*PDstandardNth23gt13 + 
        dJ223L*PDstandardNth2gt13 + J32L*J33L*PDstandardNth33gt13 + 
        dJ323L*PDstandardNth3gt13;
      
      JacPDstandardNth32gt22 = J12L*J13L*PDstandardNth11gt22 + 
        J13L*J22L*PDstandardNth12gt22 + J12L*J23L*PDstandardNth12gt22 + 
        J13L*J32L*PDstandardNth13gt22 + J12L*J33L*PDstandardNth13gt22 + 
        dJ123L*PDstandardNth1gt22 + J22L*J23L*PDstandardNth22gt22 + 
        J23L*J32L*PDstandardNth23gt22 + J22L*J33L*PDstandardNth23gt22 + 
        dJ223L*PDstandardNth2gt22 + J32L*J33L*PDstandardNth33gt22 + 
        dJ323L*PDstandardNth3gt22;
      
      JacPDstandardNth32gt23 = J12L*J13L*PDstandardNth11gt23 + 
        J13L*J22L*PDstandardNth12gt23 + J12L*J23L*PDstandardNth12gt23 + 
        J13L*J32L*PDstandardNth13gt23 + J12L*J33L*PDstandardNth13gt23 + 
        dJ123L*PDstandardNth1gt23 + J22L*J23L*PDstandardNth22gt23 + 
        J23L*J32L*PDstandardNth23gt23 + J22L*J33L*PDstandardNth23gt23 + 
        dJ223L*PDstandardNth2gt23 + J32L*J33L*PDstandardNth33gt23 + 
        dJ323L*PDstandardNth3gt23;
      
      JacPDstandardNth32gt33 = J12L*J13L*PDstandardNth11gt33 + 
        J13L*J22L*PDstandardNth12gt33 + J12L*J23L*PDstandardNth12gt33 + 
        J13L*J32L*PDstandardNth13gt33 + J12L*J33L*PDstandardNth13gt33 + 
        dJ123L*PDstandardNth1gt33 + J22L*J23L*PDstandardNth22gt33 + 
        J23L*J32L*PDstandardNth23gt33 + J22L*J33L*PDstandardNth23gt33 + 
        dJ223L*PDstandardNth2gt33 + J32L*J33L*PDstandardNth33gt33 + 
        dJ323L*PDstandardNth3gt33;
    }
    else
    {
      JacPDstandardNth1gt11 = PDstandardNth1gt11;
      
      JacPDstandardNth1gt12 = PDstandardNth1gt12;
      
      JacPDstandardNth1gt13 = PDstandardNth1gt13;
      
      JacPDstandardNth1gt22 = PDstandardNth1gt22;
      
      JacPDstandardNth1gt23 = PDstandardNth1gt23;
      
      JacPDstandardNth1gt33 = PDstandardNth1gt33;
      
      JacPDstandardNth1phi = PDstandardNth1phi;
      
      JacPDstandardNth1Xt1 = PDstandardNth1Xt1;
      
      JacPDstandardNth1Xt2 = PDstandardNth1Xt2;
      
      JacPDstandardNth1Xt3 = PDstandardNth1Xt3;
      
      JacPDstandardNth2gt11 = PDstandardNth2gt11;
      
      JacPDstandardNth2gt12 = PDstandardNth2gt12;
      
      JacPDstandardNth2gt13 = PDstandardNth2gt13;
      
      JacPDstandardNth2gt22 = PDstandardNth2gt22;
      
      JacPDstandardNth2gt23 = PDstandardNth2gt23;
      
      JacPDstandardNth2gt33 = PDstandardNth2gt33;
      
      JacPDstandardNth2phi = PDstandardNth2phi;
      
      JacPDstandardNth2Xt1 = PDstandardNth2Xt1;
      
      JacPDstandardNth2Xt2 = PDstandardNth2Xt2;
      
      JacPDstandardNth2Xt3 = PDstandardNth2Xt3;
      
      JacPDstandardNth3gt11 = PDstandardNth3gt11;
      
      JacPDstandardNth3gt12 = PDstandardNth3gt12;
      
      JacPDstandardNth3gt13 = PDstandardNth3gt13;
      
      JacPDstandardNth3gt22 = PDstandardNth3gt22;
      
      JacPDstandardNth3gt23 = PDstandardNth3gt23;
      
      JacPDstandardNth3gt33 = PDstandardNth3gt33;
      
      JacPDstandardNth3phi = PDstandardNth3phi;
      
      JacPDstandardNth3Xt1 = PDstandardNth3Xt1;
      
      JacPDstandardNth3Xt2 = PDstandardNth3Xt2;
      
      JacPDstandardNth3Xt3 = PDstandardNth3Xt3;
      
      JacPDstandardNth11gt11 = PDstandardNth11gt11;
      
      JacPDstandardNth11gt12 = PDstandardNth11gt12;
      
      JacPDstandardNth11gt13 = PDstandardNth11gt13;
      
      JacPDstandardNth11gt22 = PDstandardNth11gt22;
      
      JacPDstandardNth11gt23 = PDstandardNth11gt23;
      
      JacPDstandardNth11gt33 = PDstandardNth11gt33;
      
      JacPDstandardNth11phi = PDstandardNth11phi;
      
      JacPDstandardNth22gt11 = PDstandardNth22gt11;
      
      JacPDstandardNth22gt12 = PDstandardNth22gt12;
      
      JacPDstandardNth22gt13 = PDstandardNth22gt13;
      
      JacPDstandardNth22gt22 = PDstandardNth22gt22;
      
      JacPDstandardNth22gt23 = PDstandardNth22gt23;
      
      JacPDstandardNth22gt33 = PDstandardNth22gt33;
      
      JacPDstandardNth22phi = PDstandardNth22phi;
      
      JacPDstandardNth33gt11 = PDstandardNth33gt11;
      
      JacPDstandardNth33gt12 = PDstandardNth33gt12;
      
      JacPDstandardNth33gt13 = PDstandardNth33gt13;
      
      JacPDstandardNth33gt22 = PDstandardNth33gt22;
      
      JacPDstandardNth33gt23 = PDstandardNth33gt23;
      
      JacPDstandardNth33gt33 = PDstandardNth33gt33;
      
      JacPDstandardNth33phi = PDstandardNth33phi;
      
      JacPDstandardNth12gt11 = PDstandardNth12gt11;
      
      JacPDstandardNth12gt12 = PDstandardNth12gt12;
      
      JacPDstandardNth12gt13 = PDstandardNth12gt13;
      
      JacPDstandardNth12gt22 = PDstandardNth12gt22;
      
      JacPDstandardNth12gt23 = PDstandardNth12gt23;
      
      JacPDstandardNth12gt33 = PDstandardNth12gt33;
      
      JacPDstandardNth12phi = PDstandardNth12phi;
      
      JacPDstandardNth13gt11 = PDstandardNth13gt11;
      
      JacPDstandardNth13gt12 = PDstandardNth13gt12;
      
      JacPDstandardNth13gt13 = PDstandardNth13gt13;
      
      JacPDstandardNth13gt22 = PDstandardNth13gt22;
      
      JacPDstandardNth13gt23 = PDstandardNth13gt23;
      
      JacPDstandardNth13gt33 = PDstandardNth13gt33;
      
      JacPDstandardNth13phi = PDstandardNth13phi;
      
      JacPDstandardNth21gt11 = PDstandardNth12gt11;
      
      JacPDstandardNth21gt12 = PDstandardNth12gt12;
      
      JacPDstandardNth21gt13 = PDstandardNth12gt13;
      
      JacPDstandardNth21gt22 = PDstandardNth12gt22;
      
      JacPDstandardNth21gt23 = PDstandardNth12gt23;
      
      JacPDstandardNth21gt33 = PDstandardNth12gt33;
      
      JacPDstandardNth23gt11 = PDstandardNth23gt11;
      
      JacPDstandardNth23gt12 = PDstandardNth23gt12;
      
      JacPDstandardNth23gt13 = PDstandardNth23gt13;
      
      JacPDstandardNth23gt22 = PDstandardNth23gt22;
      
      JacPDstandardNth23gt23 = PDstandardNth23gt23;
      
      JacPDstandardNth23gt33 = PDstandardNth23gt33;
      
      JacPDstandardNth23phi = PDstandardNth23phi;
      
      JacPDstandardNth31gt11 = PDstandardNth13gt11;
      
      JacPDstandardNth31gt12 = PDstandardNth13gt12;
      
      JacPDstandardNth31gt13 = PDstandardNth13gt13;
      
      JacPDstandardNth31gt22 = PDstandardNth13gt22;
      
      JacPDstandardNth31gt23 = PDstandardNth13gt23;
      
      JacPDstandardNth31gt33 = PDstandardNth13gt33;
      
      JacPDstandardNth32gt11 = PDstandardNth23gt11;
      
      JacPDstandardNth32gt12 = PDstandardNth23gt12;
      
      JacPDstandardNth32gt13 = PDstandardNth23gt13;
      
      JacPDstandardNth32gt22 = PDstandardNth23gt22;
      
      JacPDstandardNth32gt23 = PDstandardNth23gt23;
      
      JacPDstandardNth32gt33 = PDstandardNth23gt33;
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
    
    CCTK_REAL Gtlu111 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu11 + Gtl112*gtu12 
      + Gtl113*gtu13;
    
    CCTK_REAL Gtlu112 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu12 + Gtl112*gtu22 
      + Gtl113*gtu23;
    
    CCTK_REAL Gtlu113 CCTK_ATTRIBUTE_UNUSED = Gtl111*gtu13 + Gtl112*gtu23 
      + Gtl113*gtu33;
    
    CCTK_REAL Gtlu121 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu11 + Gtl122*gtu12 
      + Gtl123*gtu13;
    
    CCTK_REAL Gtlu122 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu12 + Gtl122*gtu22 
      + Gtl123*gtu23;
    
    CCTK_REAL Gtlu123 CCTK_ATTRIBUTE_UNUSED = Gtl112*gtu13 + Gtl122*gtu23 
      + Gtl123*gtu33;
    
    CCTK_REAL Gtlu131 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu11 + Gtl123*gtu12 
      + Gtl133*gtu13;
    
    CCTK_REAL Gtlu132 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu12 + Gtl123*gtu22 
      + Gtl133*gtu23;
    
    CCTK_REAL Gtlu133 CCTK_ATTRIBUTE_UNUSED = Gtl113*gtu13 + Gtl123*gtu23 
      + Gtl133*gtu33;
    
    CCTK_REAL Gtlu211 CCTK_ATTRIBUTE_UNUSED = Gtl211*gtu11 + Gtl212*gtu12 
      + Gtl213*gtu13;
    
    CCTK_REAL Gtlu212 CCTK_ATTRIBUTE_UNUSED = Gtl211*gtu12 + Gtl212*gtu22 
      + Gtl213*gtu23;
    
    CCTK_REAL Gtlu213 CCTK_ATTRIBUTE_UNUSED = Gtl211*gtu13 + Gtl212*gtu23 
      + Gtl213*gtu33;
    
    CCTK_REAL Gtlu221 CCTK_ATTRIBUTE_UNUSED = Gtl212*gtu11 + Gtl222*gtu12 
      + Gtl223*gtu13;
    
    CCTK_REAL Gtlu222 CCTK_ATTRIBUTE_UNUSED = Gtl212*gtu12 + Gtl222*gtu22 
      + Gtl223*gtu23;
    
    CCTK_REAL Gtlu223 CCTK_ATTRIBUTE_UNUSED = Gtl212*gtu13 + Gtl222*gtu23 
      + Gtl223*gtu33;
    
    CCTK_REAL Gtlu231 CCTK_ATTRIBUTE_UNUSED = Gtl213*gtu11 + Gtl223*gtu12 
      + Gtl233*gtu13;
    
    CCTK_REAL Gtlu232 CCTK_ATTRIBUTE_UNUSED = Gtl213*gtu12 + Gtl223*gtu22 
      + Gtl233*gtu23;
    
    CCTK_REAL Gtlu233 CCTK_ATTRIBUTE_UNUSED = Gtl213*gtu13 + Gtl223*gtu23 
      + Gtl233*gtu33;
    
    CCTK_REAL Gtlu311 CCTK_ATTRIBUTE_UNUSED = Gtl311*gtu11 + Gtl312*gtu12 
      + Gtl313*gtu13;
    
    CCTK_REAL Gtlu312 CCTK_ATTRIBUTE_UNUSED = Gtl311*gtu12 + Gtl312*gtu22 
      + Gtl313*gtu23;
    
    CCTK_REAL Gtlu313 CCTK_ATTRIBUTE_UNUSED = Gtl311*gtu13 + Gtl312*gtu23 
      + Gtl313*gtu33;
    
    CCTK_REAL Gtlu321 CCTK_ATTRIBUTE_UNUSED = Gtl312*gtu11 + Gtl322*gtu12 
      + Gtl323*gtu13;
    
    CCTK_REAL Gtlu322 CCTK_ATTRIBUTE_UNUSED = Gtl312*gtu12 + Gtl322*gtu22 
      + Gtl323*gtu23;
    
    CCTK_REAL Gtlu323 CCTK_ATTRIBUTE_UNUSED = Gtl312*gtu13 + Gtl322*gtu23 
      + Gtl323*gtu33;
    
    CCTK_REAL Gtlu331 CCTK_ATTRIBUTE_UNUSED = Gtl313*gtu11 + Gtl323*gtu12 
      + Gtl333*gtu13;
    
    CCTK_REAL Gtlu332 CCTK_ATTRIBUTE_UNUSED = Gtl313*gtu12 + Gtl323*gtu22 
      + Gtl333*gtu23;
    
    CCTK_REAL Gtlu333 CCTK_ATTRIBUTE_UNUSED = Gtl313*gtu13 + Gtl323*gtu23 
      + Gtl333*gtu33;
    
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
    
    CCTK_REAL gu11 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu11;
    
    CCTK_REAL gu12 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu12;
    
    CCTK_REAL gu13 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu13;
    
    CCTK_REAL gu22 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu22;
    
    CCTK_REAL gu23 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu23;
    
    CCTK_REAL gu33 CCTK_ATTRIBUTE_UNUSED = em4phi*gtu33;
    
    CCTK_REAL Rt11 CCTK_ATTRIBUTE_UNUSED = 3*Gt111*Gtlu111 + 
      3*Gt112*Gtlu112 + 3*Gt113*Gtlu113 + 2*Gt211*Gtlu121 + 2*Gt212*Gtlu122 + 
      2*Gt213*Gtlu123 + 2*Gt311*Gtlu131 + 2*Gt312*Gtlu132 + 2*Gt313*Gtlu133 + 
      Gt211*Gtlu211 + Gt212*Gtlu212 + Gt213*Gtlu213 + Gt311*Gtlu311 + 
      Gt312*Gtlu312 + Gt313*Gtlu313 + gt11L*JacPDstandardNth1Xt1 + 
      gt12L*JacPDstandardNth1Xt2 + gt13L*JacPDstandardNth1Xt3 + 
      0.5*(-(gtu11*JacPDstandardNth11gt11) - gtu12*JacPDstandardNth12gt11 - 
      gtu13*JacPDstandardNth13gt11 - gtu12*JacPDstandardNth21gt11 - 
      gtu22*JacPDstandardNth22gt11 - gtu23*JacPDstandardNth23gt11 - 
      gtu13*JacPDstandardNth31gt11 - gtu23*JacPDstandardNth32gt11 - 
      gtu33*JacPDstandardNth33gt11) + Gtl111*Xtn1 + Gtl112*Xtn2 + 
      Gtl113*Xtn3;
    
    CCTK_REAL Rt12 CCTK_ATTRIBUTE_UNUSED = Gt112*Gtlu111 + Gt122*Gtlu112 + 
      Gt123*Gtlu113 + Gt111*Gtlu121 + Gt212*Gtlu121 + Gt112*Gtlu122 + 
      Gt222*Gtlu122 + Gt113*Gtlu123 + Gt223*Gtlu123 + Gt312*Gtlu131 + 
      Gt322*Gtlu132 + Gt323*Gtlu133 + Gt111*Gtlu211 + Gt112*Gtlu212 + 
      Gt113*Gtlu213 + 2*Gt211*Gtlu221 + 2*Gt212*Gtlu222 + 2*Gt213*Gtlu223 + 
      Gt311*Gtlu231 + Gt312*Gtlu232 + Gt313*Gtlu233 + Gt311*Gtlu321 + 
      Gt312*Gtlu322 + Gt313*Gtlu323 + 0.5*(gt12L*JacPDstandardNth1Xt1 + 
      gt22L*JacPDstandardNth1Xt2 + gt23L*JacPDstandardNth1Xt3) + 
      0.5*(gt11L*JacPDstandardNth2Xt1 + gt12L*JacPDstandardNth2Xt2 + 
      gt13L*JacPDstandardNth2Xt3) + 0.5*(-(gtu11*JacPDstandardNth11gt12) - 
      gtu12*JacPDstandardNth12gt12 - gtu13*JacPDstandardNth13gt12 - 
      gtu12*JacPDstandardNth21gt12 - gtu22*JacPDstandardNth22gt12 - 
      gtu23*JacPDstandardNth23gt12 - gtu13*JacPDstandardNth31gt12 - 
      gtu23*JacPDstandardNth32gt12 - gtu33*JacPDstandardNth33gt12) + 
      0.5*(Gtl112*Xtn1 + Gtl122*Xtn2 + Gtl123*Xtn3) + 0.5*(Gtl211*Xtn1 + 
      Gtl212*Xtn2 + Gtl213*Xtn3);
    
    CCTK_REAL Rt13 CCTK_ATTRIBUTE_UNUSED = Gt113*Gtlu111 + Gt123*Gtlu112 + 
      Gt133*Gtlu113 + Gt213*Gtlu121 + Gt223*Gtlu122 + Gt233*Gtlu123 + 
      Gt111*Gtlu131 + Gt313*Gtlu131 + Gt112*Gtlu132 + Gt323*Gtlu132 + 
      Gt113*Gtlu133 + Gt333*Gtlu133 + Gt211*Gtlu231 + Gt212*Gtlu232 + 
      Gt213*Gtlu233 + Gt111*Gtlu311 + Gt112*Gtlu312 + Gt113*Gtlu313 + 
      Gt211*Gtlu321 + Gt212*Gtlu322 + Gt213*Gtlu323 + 2*Gt311*Gtlu331 + 
      2*Gt312*Gtlu332 + 2*Gt313*Gtlu333 + 0.5*(gt13L*JacPDstandardNth1Xt1 + 
      gt23L*JacPDstandardNth1Xt2 + gt33L*JacPDstandardNth1Xt3) + 
      0.5*(-(gtu11*JacPDstandardNth11gt13) - gtu12*JacPDstandardNth12gt13 - 
      gtu13*JacPDstandardNth13gt13 - gtu12*JacPDstandardNth21gt13 - 
      gtu22*JacPDstandardNth22gt13 - gtu23*JacPDstandardNth23gt13 - 
      gtu13*JacPDstandardNth31gt13 - gtu23*JacPDstandardNth32gt13 - 
      gtu33*JacPDstandardNth33gt13) + 0.5*(gt11L*JacPDstandardNth3Xt1 + 
      gt12L*JacPDstandardNth3Xt2 + gt13L*JacPDstandardNth3Xt3) + 
      0.5*(Gtl113*Xtn1 + Gtl123*Xtn2 + Gtl133*Xtn3) + 0.5*(Gtl311*Xtn1 + 
      Gtl312*Xtn2 + Gtl313*Xtn3);
    
    CCTK_REAL Rt22 CCTK_ATTRIBUTE_UNUSED = Gt112*Gtlu121 + Gt122*Gtlu122 + 
      Gt123*Gtlu123 + 2*Gt112*Gtlu211 + 2*Gt122*Gtlu212 + 2*Gt123*Gtlu213 + 
      3*Gt212*Gtlu221 + 3*Gt222*Gtlu222 + 3*Gt223*Gtlu223 + 2*Gt312*Gtlu231 + 
      2*Gt322*Gtlu232 + 2*Gt323*Gtlu233 + Gt312*Gtlu321 + Gt322*Gtlu322 + 
      Gt323*Gtlu323 + gt12L*JacPDstandardNth2Xt1 + gt22L*JacPDstandardNth2Xt2 
      + gt23L*JacPDstandardNth2Xt3 + 0.5*(-(gtu11*JacPDstandardNth11gt22) - 
      gtu12*JacPDstandardNth12gt22 - gtu13*JacPDstandardNth13gt22 - 
      gtu12*JacPDstandardNth21gt22 - gtu22*JacPDstandardNth22gt22 - 
      gtu23*JacPDstandardNth23gt22 - gtu13*JacPDstandardNth31gt22 - 
      gtu23*JacPDstandardNth32gt22 - gtu33*JacPDstandardNth33gt22) + 
      Gtl212*Xtn1 + Gtl222*Xtn2 + Gtl223*Xtn3;
    
    CCTK_REAL Rt23 CCTK_ATTRIBUTE_UNUSED = Gt112*Gtlu131 + Gt122*Gtlu132 + 
      Gt123*Gtlu133 + Gt113*Gtlu211 + Gt123*Gtlu212 + Gt133*Gtlu213 + 
      Gt213*Gtlu221 + Gt223*Gtlu222 + Gt233*Gtlu223 + Gt212*Gtlu231 + 
      Gt313*Gtlu231 + Gt222*Gtlu232 + Gt323*Gtlu232 + Gt223*Gtlu233 + 
      Gt333*Gtlu233 + Gt112*Gtlu311 + Gt122*Gtlu312 + Gt123*Gtlu313 + 
      Gt212*Gtlu321 + Gt222*Gtlu322 + Gt223*Gtlu323 + 2*Gt312*Gtlu331 + 
      2*Gt322*Gtlu332 + 2*Gt323*Gtlu333 + 0.5*(gt13L*JacPDstandardNth2Xt1 + 
      gt23L*JacPDstandardNth2Xt2 + gt33L*JacPDstandardNth2Xt3) + 
      0.5*(-(gtu11*JacPDstandardNth11gt23) - gtu12*JacPDstandardNth12gt23 - 
      gtu13*JacPDstandardNth13gt23 - gtu12*JacPDstandardNth21gt23 - 
      gtu22*JacPDstandardNth22gt23 - gtu23*JacPDstandardNth23gt23 - 
      gtu13*JacPDstandardNth31gt23 - gtu23*JacPDstandardNth32gt23 - 
      gtu33*JacPDstandardNth33gt23) + 0.5*(gt12L*JacPDstandardNth3Xt1 + 
      gt22L*JacPDstandardNth3Xt2 + gt23L*JacPDstandardNth3Xt3) + 
      0.5*(Gtl213*Xtn1 + Gtl223*Xtn2 + Gtl233*Xtn3) + 0.5*(Gtl312*Xtn1 + 
      Gtl322*Xtn2 + Gtl323*Xtn3);
    
    CCTK_REAL Rt33 CCTK_ATTRIBUTE_UNUSED = Gt113*Gtlu131 + Gt123*Gtlu132 + 
      Gt133*Gtlu133 + Gt213*Gtlu231 + Gt223*Gtlu232 + Gt233*Gtlu233 + 
      2*Gt113*Gtlu311 + 2*Gt123*Gtlu312 + 2*Gt133*Gtlu313 + 2*Gt213*Gtlu321 + 
      2*Gt223*Gtlu322 + 2*Gt233*Gtlu323 + 3*Gt313*Gtlu331 + 3*Gt323*Gtlu332 + 
      3*Gt333*Gtlu333 + 0.5*(-(gtu11*JacPDstandardNth11gt33) - 
      gtu12*JacPDstandardNth12gt33 - gtu13*JacPDstandardNth13gt33 - 
      gtu12*JacPDstandardNth21gt33 - gtu22*JacPDstandardNth22gt33 - 
      gtu23*JacPDstandardNth23gt33 - gtu13*JacPDstandardNth31gt33 - 
      gtu23*JacPDstandardNth32gt33 - gtu33*JacPDstandardNth33gt33) + 
      gt13L*JacPDstandardNth3Xt1 + gt23L*JacPDstandardNth3Xt2 + 
      gt33L*JacPDstandardNth3Xt3 + Gtl313*Xtn1 + Gtl323*Xtn2 + Gtl333*Xtn3;
    
    CCTK_REAL fac1 CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod == 
      1,-0.5*INV(phiL),1);
    
    CCTK_REAL cdphi1 CCTK_ATTRIBUTE_UNUSED = fac1*JacPDstandardNth1phi;
    
    CCTK_REAL cdphi2 CCTK_ATTRIBUTE_UNUSED = fac1*JacPDstandardNth2phi;
    
    CCTK_REAL cdphi3 CCTK_ATTRIBUTE_UNUSED = fac1*JacPDstandardNth3phi;
    
    CCTK_REAL fac2 CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod == 
      1,0.5*INV(SQR(phiL)),0);
    
    CCTK_REAL cdphi211 CCTK_ATTRIBUTE_UNUSED = fac1*(JacPDstandardNth11phi 
      - Gt111*JacPDstandardNth1phi - Gt211*JacPDstandardNth2phi - 
      Gt311*JacPDstandardNth3phi) + fac2*SQR(JacPDstandardNth1phi);
    
    CCTK_REAL cdphi212 CCTK_ATTRIBUTE_UNUSED = 
      fac2*JacPDstandardNth1phi*JacPDstandardNth2phi + 
      fac1*(JacPDstandardNth12phi - Gt112*JacPDstandardNth1phi - 
      Gt212*JacPDstandardNth2phi - Gt312*JacPDstandardNth3phi);
    
    CCTK_REAL cdphi213 CCTK_ATTRIBUTE_UNUSED = 
      fac2*JacPDstandardNth1phi*JacPDstandardNth3phi + 
      fac1*(JacPDstandardNth13phi - Gt113*JacPDstandardNth1phi - 
      Gt213*JacPDstandardNth2phi - Gt313*JacPDstandardNth3phi);
    
    CCTK_REAL cdphi222 CCTK_ATTRIBUTE_UNUSED = 
      fac1*(-(Gt122*JacPDstandardNth1phi) + JacPDstandardNth22phi - 
      Gt222*JacPDstandardNth2phi - Gt322*JacPDstandardNth3phi) + 
      fac2*SQR(JacPDstandardNth2phi);
    
    CCTK_REAL cdphi223 CCTK_ATTRIBUTE_UNUSED = 
      fac2*JacPDstandardNth2phi*JacPDstandardNth3phi + 
      fac1*(-(Gt123*JacPDstandardNth1phi) + JacPDstandardNth23phi - 
      Gt223*JacPDstandardNth2phi - Gt323*JacPDstandardNth3phi);
    
    CCTK_REAL cdphi233 CCTK_ATTRIBUTE_UNUSED = 
      fac1*(-(Gt133*JacPDstandardNth1phi) - Gt233*JacPDstandardNth2phi + 
      JacPDstandardNth33phi - Gt333*JacPDstandardNth3phi) + 
      fac2*SQR(JacPDstandardNth3phi);
    
    CCTK_REAL Rphi11 CCTK_ATTRIBUTE_UNUSED = -2*cdphi211 - 
      2*gt11L*(cdphi211*gtu11 + 2*cdphi212*gtu12 + 2*cdphi213*gtu13 + 
      cdphi222*gtu22 + 2*cdphi223*gtu23 + cdphi233*gtu33) - 
      4*gt11L*(cdphi1*(cdphi1*gtu11 + cdphi2*gtu12 + cdphi3*gtu13) + 
      cdphi2*(cdphi1*gtu12 + cdphi2*gtu22 + cdphi3*gtu23) + 
      cdphi3*(cdphi1*gtu13 + cdphi2*gtu23 + cdphi3*gtu33)) + 4*SQR(cdphi1);
    
    CCTK_REAL Rphi12 CCTK_ATTRIBUTE_UNUSED = 4*cdphi1*cdphi2 - 2*cdphi212 
      - 2*gt12L*(cdphi211*gtu11 + 2*cdphi212*gtu12 + 2*cdphi213*gtu13 + 
      cdphi222*gtu22 + 2*cdphi223*gtu23 + cdphi233*gtu33) - 
      4*gt12L*(cdphi1*(cdphi1*gtu11 + cdphi2*gtu12 + cdphi3*gtu13) + 
      cdphi2*(cdphi1*gtu12 + cdphi2*gtu22 + cdphi3*gtu23) + 
      cdphi3*(cdphi1*gtu13 + cdphi2*gtu23 + cdphi3*gtu33));
    
    CCTK_REAL Rphi13 CCTK_ATTRIBUTE_UNUSED = -2*cdphi213 + 4*cdphi1*cdphi3 
      - 2*gt13L*(cdphi211*gtu11 + 2*cdphi212*gtu12 + 2*cdphi213*gtu13 + 
      cdphi222*gtu22 + 2*cdphi223*gtu23 + cdphi233*gtu33) - 
      4*gt13L*(cdphi1*(cdphi1*gtu11 + cdphi2*gtu12 + cdphi3*gtu13) + 
      cdphi2*(cdphi1*gtu12 + cdphi2*gtu22 + cdphi3*gtu23) + 
      cdphi3*(cdphi1*gtu13 + cdphi2*gtu23 + cdphi3*gtu33));
    
    CCTK_REAL Rphi22 CCTK_ATTRIBUTE_UNUSED = -2*cdphi222 - 
      2*gt22L*(cdphi211*gtu11 + 2*cdphi212*gtu12 + 2*cdphi213*gtu13 + 
      cdphi222*gtu22 + 2*cdphi223*gtu23 + cdphi233*gtu33) - 
      4*gt22L*(cdphi1*(cdphi1*gtu11 + cdphi2*gtu12 + cdphi3*gtu13) + 
      cdphi2*(cdphi1*gtu12 + cdphi2*gtu22 + cdphi3*gtu23) + 
      cdphi3*(cdphi1*gtu13 + cdphi2*gtu23 + cdphi3*gtu33)) + 4*SQR(cdphi2);
    
    CCTK_REAL Rphi23 CCTK_ATTRIBUTE_UNUSED = -2*cdphi223 + 4*cdphi2*cdphi3 
      - 2*gt23L*(cdphi211*gtu11 + 2*cdphi212*gtu12 + 2*cdphi213*gtu13 + 
      cdphi222*gtu22 + 2*cdphi223*gtu23 + cdphi233*gtu33) - 
      4*gt23L*(cdphi1*(cdphi1*gtu11 + cdphi2*gtu12 + cdphi3*gtu13) + 
      cdphi2*(cdphi1*gtu12 + cdphi2*gtu22 + cdphi3*gtu23) + 
      cdphi3*(cdphi1*gtu13 + cdphi2*gtu23 + cdphi3*gtu33));
    
    CCTK_REAL Rphi33 CCTK_ATTRIBUTE_UNUSED = -2*cdphi233 - 
      2*gt33L*(cdphi211*gtu11 + 2*cdphi212*gtu12 + 2*cdphi213*gtu13 + 
      cdphi222*gtu22 + 2*cdphi223*gtu23 + cdphi233*gtu33) - 
      4*gt33L*(cdphi1*(cdphi1*gtu11 + cdphi2*gtu12 + cdphi3*gtu13) + 
      cdphi2*(cdphi1*gtu12 + cdphi2*gtu22 + cdphi3*gtu23) + 
      cdphi3*(cdphi1*gtu13 + cdphi2*gtu23 + cdphi3*gtu33)) + 4*SQR(cdphi3);
    
    CCTK_REAL R11 CCTK_ATTRIBUTE_UNUSED = Rphi11 + Rt11;
    
    CCTK_REAL R12 CCTK_ATTRIBUTE_UNUSED = Rphi12 + Rt12;
    
    CCTK_REAL R13 CCTK_ATTRIBUTE_UNUSED = Rphi13 + Rt13;
    
    CCTK_REAL R22 CCTK_ATTRIBUTE_UNUSED = Rphi22 + Rt22;
    
    CCTK_REAL R23 CCTK_ATTRIBUTE_UNUSED = Rphi23 + Rt23;
    
    CCTK_REAL R33 CCTK_ATTRIBUTE_UNUSED = Rphi33 + Rt33;
    
    CCTK_REAL trR CCTK_ATTRIBUTE_UNUSED = gu11*R11 + 2*gu12*R12 + 
      2*gu13*R13 + gu22*R22 + 2*gu23*R23 + gu33*R33;
    
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
    
    CCTK_REAL rho CCTK_ATTRIBUTE_UNUSED = (eTttL - 2*(beta1L*eTtxL + 
      beta2L*eTtyL + beta3L*eTtzL) + beta1L*(beta1L*eTxxL + beta2L*eTxyL + 
      beta3L*eTxzL) + beta2L*(beta1L*eTxyL + beta2L*eTyyL + beta3L*eTyzL) + 
      beta3L*(beta1L*eTxzL + beta2L*eTyzL + beta3L*eTzzL))*INV(SQR(alphaL));
    
    CCTK_REAL HL CCTK_ATTRIBUTE_UNUSED = -2*Atm12*Atm21 - 2*Atm13*Atm31 - 
      2*Atm23*Atm32 - 16*Pi*rho + trR + 
      0.666666666666666666666666666667*SQR(trKL) - SQR(Atm11) - SQR(Atm22) - 
      SQR(Atm33);
    /* Copy local copies back to grid functions */
    H[index] = HL;
  }
  CCTK_ENDLOOP3(ML_BSSN_NoVec_constraints1);
}
extern "C" void ML_BSSN_NoVec_constraints1(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_NoVec_constraints1_Body");
  }
  if (cctk_iteration % ML_BSSN_NoVec_constraints1_calc_every != ML_BSSN_NoVec_constraints1_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ML_BSSN_NoVec::ML_curv",
    "ML_BSSN_NoVec::ML_Gamma",
    "ML_BSSN_NoVec::ML_Ham",
    "ML_BSSN_NoVec::ML_lapse",
    "ML_BSSN_NoVec::ML_log_confac",
    "ML_BSSN_NoVec::ML_metric",
    "ML_BSSN_NoVec::ML_shift",
    "ML_BSSN_NoVec::ML_trace_curv"};
  AssertGroupStorage(cctkGH, "ML_BSSN_NoVec_constraints1", 8, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_NoVec_constraints1", 1, 1, 1);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_NoVec_constraints1", 2, 2, 2);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, ML_BSSN_NoVec_constraints1_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_NoVec_constraints1_Body");
  }
}

} // namespace ML_BSSN_NoVec
