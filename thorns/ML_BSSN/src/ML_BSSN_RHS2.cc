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

namespace ML_BSSN {

extern "C" void ML_BSSN_RHS2_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_RHS2_calc_every != ML_BSSN_RHS2_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_curvrhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN::ML_curvrhs.");
  return;
}

static void ML_BSSN_RHS2_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Include user-supplied include files */
  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const ptrdiff_t cctkLbnd1 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[0];
  const ptrdiff_t cctkLbnd2 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[1];
  const ptrdiff_t cctkLbnd3 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[2];
  const CCTK_REAL_VEC t CCTK_ATTRIBUTE_UNUSED = ToReal(cctk_time);
  const CCTK_REAL_VEC cctkOriginSpace1 CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_ORIGIN_SPACE(0));
  const CCTK_REAL_VEC cctkOriginSpace2 CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_ORIGIN_SPACE(1));
  const CCTK_REAL_VEC cctkOriginSpace3 CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_ORIGIN_SPACE(2));
  const CCTK_REAL_VEC dt CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_TIME);
  const CCTK_REAL_VEC dx CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_SPACE(0));
  const CCTK_REAL_VEC dy CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_SPACE(1));
  const CCTK_REAL_VEC dz CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_SPACE(2));
  const CCTK_REAL_VEC dxi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dx);
  const CCTK_REAL_VEC dyi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dy);
  const CCTK_REAL_VEC dzi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dz);
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
  const CCTK_REAL_VEC p1o144dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dx,dy));
  const CCTK_REAL_VEC p1o144dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dx,dz));
  const CCTK_REAL_VEC p1o144dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dy,dz));
  const CCTK_REAL_VEC p1o24dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dx);
  const CCTK_REAL_VEC p1o24dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dy);
  const CCTK_REAL_VEC p1o24dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dz);
  const CCTK_REAL_VEC p1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dx);
  const CCTK_REAL_VEC p1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dy);
  const CCTK_REAL_VEC p1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dz);
  const CCTK_REAL_VEC p1o4dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dx);
  const CCTK_REAL_VEC p1o4dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dx,dy));
  const CCTK_REAL_VEC p1o4dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dx,dz));
  const CCTK_REAL_VEC p1o4dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dy);
  const CCTK_REAL_VEC p1o4dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dy,dz));
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
  const bool usejacobian1 = (!CCTK_IsFunctionAliased("MultiPatch_GetMap") || MultiPatch_GetMap(cctkGH) != jacobian_identity_map)
                        && strlen(jacobian_group) > 0;
  const bool usejacobian = assume_use_jacobian>=0 ? assume_use_jacobian : usejacobian1;
  if (usejacobian && (strlen(jacobian_derivative_group) == 0))
  {
    CCTK_WARN(1, "GenericFD::jacobian_group and GenericFD::jacobian_derivative_group must both be set to valid group names");
  }
  
  const CCTK_REAL* restrict jacobian_ptrs[9];
  if (usejacobian) GroupDataPointers(cctkGH, jacobian_group,
                                                9, jacobian_ptrs);
  
  const CCTK_REAL* restrict const J11 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[0] : 0;
  const CCTK_REAL* restrict const J12 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[1] : 0;
  const CCTK_REAL* restrict const J13 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[2] : 0;
  const CCTK_REAL* restrict const J21 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[3] : 0;
  const CCTK_REAL* restrict const J22 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[4] : 0;
  const CCTK_REAL* restrict const J23 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[5] : 0;
  const CCTK_REAL* restrict const J31 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[6] : 0;
  const CCTK_REAL* restrict const J32 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[7] : 0;
  const CCTK_REAL* restrict const J33 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[8] : 0;
  
  const CCTK_REAL* restrict jacobian_determinant_ptrs[1] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian && strlen(jacobian_determinant_group) > 0) GroupDataPointers(cctkGH, jacobian_determinant_group,
                                                1, jacobian_determinant_ptrs);
  
  const CCTK_REAL* restrict const detJ CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_determinant_ptrs[0] : 0;
  
  const CCTK_REAL* restrict jacobian_inverse_ptrs[9] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian && strlen(jacobian_inverse_group) > 0) GroupDataPointers(cctkGH, jacobian_inverse_group,
                                                9, jacobian_inverse_ptrs);
  
  const CCTK_REAL* restrict const iJ11 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[0] : 0;
  const CCTK_REAL* restrict const iJ12 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[1] : 0;
  const CCTK_REAL* restrict const iJ13 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[2] : 0;
  const CCTK_REAL* restrict const iJ21 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[3] : 0;
  const CCTK_REAL* restrict const iJ22 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[4] : 0;
  const CCTK_REAL* restrict const iJ23 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[5] : 0;
  const CCTK_REAL* restrict const iJ31 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[6] : 0;
  const CCTK_REAL* restrict const iJ32 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[7] : 0;
  const CCTK_REAL* restrict const iJ33 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[8] : 0;
  
  const CCTK_REAL* restrict jacobian_derivative_ptrs[18] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian) GroupDataPointers(cctkGH, jacobian_derivative_group,
                                      18, jacobian_derivative_ptrs);
  
  const CCTK_REAL* restrict const dJ111 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[0] : 0;
  const CCTK_REAL* restrict const dJ112 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[1] : 0;
  const CCTK_REAL* restrict const dJ113 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[2] : 0;
  const CCTK_REAL* restrict const dJ122 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[3] : 0;
  const CCTK_REAL* restrict const dJ123 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[4] : 0;
  const CCTK_REAL* restrict const dJ133 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[5] : 0;
  const CCTK_REAL* restrict const dJ211 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[6] : 0;
  const CCTK_REAL* restrict const dJ212 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[7] : 0;
  const CCTK_REAL* restrict const dJ213 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[8] : 0;
  const CCTK_REAL* restrict const dJ222 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[9] : 0;
  const CCTK_REAL* restrict const dJ223 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[10] : 0;
  const CCTK_REAL* restrict const dJ233 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[11] : 0;
  const CCTK_REAL* restrict const dJ311 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[12] : 0;
  const CCTK_REAL* restrict const dJ312 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[13] : 0;
  const CCTK_REAL* restrict const dJ313 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[14] : 0;
  const CCTK_REAL* restrict const dJ322 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[15] : 0;
  const CCTK_REAL* restrict const dJ323 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[16] : 0;
  const CCTK_REAL* restrict const dJ333 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[17] : 0;
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
  CCTK_LOOP3STR(ML_BSSN_RHS2,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC alphaL CCTK_ATTRIBUTE_UNUSED = vec_load(alpha[index]);
    CCTK_REAL_VEC At11L CCTK_ATTRIBUTE_UNUSED = vec_load(At11[index]);
    CCTK_REAL_VEC At12L CCTK_ATTRIBUTE_UNUSED = vec_load(At12[index]);
    CCTK_REAL_VEC At13L CCTK_ATTRIBUTE_UNUSED = vec_load(At13[index]);
    CCTK_REAL_VEC At22L CCTK_ATTRIBUTE_UNUSED = vec_load(At22[index]);
    CCTK_REAL_VEC At23L CCTK_ATTRIBUTE_UNUSED = vec_load(At23[index]);
    CCTK_REAL_VEC At33L CCTK_ATTRIBUTE_UNUSED = vec_load(At33[index]);
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
    CCTK_REAL_VEC trKL CCTK_ATTRIBUTE_UNUSED = vec_load(trK[index]);
    CCTK_REAL_VEC Xt1L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt1[index]);
    CCTK_REAL_VEC Xt2L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt2[index]);
    CCTK_REAL_VEC Xt3L CCTK_ATTRIBUTE_UNUSED = vec_load(Xt3[index]);
    
    CCTK_REAL_VEC eTxxL, eTxyL, eTxzL, eTyyL, eTyzL, eTzzL CCTK_ATTRIBUTE_UNUSED ;
    
    if (assume_stress_energy_state>=0 ? assume_stress_energy_state : *stress_energy_state)
    {
      eTxxL = vec_load(eTxx[index]);
      eTxyL = vec_load(eTxy[index]);
      eTxzL = vec_load(eTxz[index]);
      eTyyL = vec_load(eTyy[index]);
      eTyzL = vec_load(eTyz[index]);
      eTzzL = vec_load(eTzz[index]);
    }
    else
    {
      eTxxL = ToReal(0.0);
      eTxyL = ToReal(0.0);
      eTxzL = ToReal(0.0);
      eTyyL = ToReal(0.0);
      eTyzL = ToReal(0.0);
      eTzzL = ToReal(0.0);
    }
    
    CCTK_REAL_VEC dJ111L, dJ112L, dJ113L, dJ122L, dJ123L, dJ133L, dJ211L, dJ212L, dJ213L, dJ222L, dJ223L, dJ233L, dJ311L, dJ312L, dJ313L, dJ322L, dJ323L, dJ333L, J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L CCTK_ATTRIBUTE_UNUSED ;
    
    if (usejacobian)
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
    CCTK_REAL_VEC PDstandardNth1beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3beta2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3beta3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth11phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth22phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth33phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth12phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth13phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth23phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC PDstandardNth3Xt3 CCTK_ATTRIBUTE_UNUSED;
    
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
        PDstandardNth1beta2 = PDstandardNthfdOrder21(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder22(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder23(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder21(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder22(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder23(&beta3[index]);
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
        PDstandardNth1beta2 = PDstandardNthfdOrder41(&beta2[index]);
        PDstandardNth2beta2 = PDstandardNthfdOrder42(&beta2[index]);
        PDstandardNth3beta2 = PDstandardNthfdOrder43(&beta2[index]);
        PDstandardNth1beta3 = PDstandardNthfdOrder41(&beta3[index]);
        PDstandardNth2beta3 = PDstandardNthfdOrder42(&beta3[index]);
        PDstandardNth3beta3 = PDstandardNthfdOrder43(&beta3[index]);
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
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = kisgn(beta1L);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = kisgn(beta2L);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = kisgn(beta3L);
    
    CCTK_REAL_VEC JacPDstandardNth11alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth11gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth11gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth11gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth11gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth11gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth11gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth11phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth12phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth13phi CCTK_ATTRIBUTE_UNUSED;
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
    CCTK_REAL_VEC JacPDstandardNth1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth21gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth22phi CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth23phi CCTK_ATTRIBUTE_UNUSED;
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
    CCTK_REAL_VEC JacPDstandardNth2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth31gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth32gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33alpha CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33gt11 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33gt12 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33gt13 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33gt22 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33gt23 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33gt33 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth33phi CCTK_ATTRIBUTE_UNUSED;
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
    CCTK_REAL_VEC JacPDstandardNth3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC JacPDstandardNth3Xt3 CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
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
      
      JacPDstandardNth1Xt1 = 
        kmadd(J11L,PDstandardNth1Xt1,kmadd(J21L,PDstandardNth2Xt1,kmul(J31L,PDstandardNth3Xt1)));
      
      JacPDstandardNth1Xt2 = 
        kmadd(J11L,PDstandardNth1Xt2,kmadd(J21L,PDstandardNth2Xt2,kmul(J31L,PDstandardNth3Xt2)));
      
      JacPDstandardNth1Xt3 = 
        kmadd(J11L,PDstandardNth1Xt3,kmadd(J21L,PDstandardNth2Xt3,kmul(J31L,PDstandardNth3Xt3)));
      
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
      
      JacPDstandardNth2Xt1 = 
        kmadd(J12L,PDstandardNth1Xt1,kmadd(J22L,PDstandardNth2Xt1,kmul(J32L,PDstandardNth3Xt1)));
      
      JacPDstandardNth2Xt2 = 
        kmadd(J12L,PDstandardNth1Xt2,kmadd(J22L,PDstandardNth2Xt2,kmul(J32L,PDstandardNth3Xt2)));
      
      JacPDstandardNth2Xt3 = 
        kmadd(J12L,PDstandardNth1Xt3,kmadd(J22L,PDstandardNth2Xt3,kmul(J32L,PDstandardNth3Xt3)));
      
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
      
      JacPDstandardNth3Xt1 = 
        kmadd(J13L,PDstandardNth1Xt1,kmadd(J23L,PDstandardNth2Xt1,kmul(J33L,PDstandardNth3Xt1)));
      
      JacPDstandardNth3Xt2 = 
        kmadd(J13L,PDstandardNth1Xt2,kmadd(J23L,PDstandardNth2Xt2,kmul(J33L,PDstandardNth3Xt2)));
      
      JacPDstandardNth3Xt3 = 
        kmadd(J13L,PDstandardNth1Xt3,kmadd(J23L,PDstandardNth2Xt3,kmul(J33L,PDstandardNth3Xt3)));
      
      JacPDstandardNth11alpha = 
        kmadd(ToReal(2),kmul(J11L,kmul(J21L,PDstandardNth12alpha)),kmadd(ToReal(2),kmul(J11L,kmul(J31L,PDstandardNth13alpha)),kmadd(dJ111L,PDstandardNth1alpha,kmadd(ToReal(2),kmul(J21L,kmul(J31L,PDstandardNth23alpha)),kmadd(dJ211L,PDstandardNth2alpha,kmadd(dJ311L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,kmul(J11L,J11L),kmadd(PDstandardNth22alpha,kmul(J21L,J21L),kmul(PDstandardNth33alpha,kmul(J31L,J31L))))))))));
      
      JacPDstandardNth11gt11 = 
        kmadd(ToReal(2),kmul(J11L,kmul(J21L,PDstandardNth12gt11)),kmadd(ToReal(2),kmul(J11L,kmul(J31L,PDstandardNth13gt11)),kmadd(dJ111L,PDstandardNth1gt11,kmadd(ToReal(2),kmul(J21L,kmul(J31L,PDstandardNth23gt11)),kmadd(dJ211L,PDstandardNth2gt11,kmadd(dJ311L,PDstandardNth3gt11,kmadd(PDstandardNth11gt11,kmul(J11L,J11L),kmadd(PDstandardNth22gt11,kmul(J21L,J21L),kmul(PDstandardNth33gt11,kmul(J31L,J31L))))))))));
      
      JacPDstandardNth11gt12 = 
        kmadd(ToReal(2),kmul(J11L,kmul(J21L,PDstandardNth12gt12)),kmadd(ToReal(2),kmul(J11L,kmul(J31L,PDstandardNth13gt12)),kmadd(dJ111L,PDstandardNth1gt12,kmadd(ToReal(2),kmul(J21L,kmul(J31L,PDstandardNth23gt12)),kmadd(dJ211L,PDstandardNth2gt12,kmadd(dJ311L,PDstandardNth3gt12,kmadd(PDstandardNth11gt12,kmul(J11L,J11L),kmadd(PDstandardNth22gt12,kmul(J21L,J21L),kmul(PDstandardNth33gt12,kmul(J31L,J31L))))))))));
      
      JacPDstandardNth11gt13 = 
        kmadd(ToReal(2),kmul(J11L,kmul(J21L,PDstandardNth12gt13)),kmadd(ToReal(2),kmul(J11L,kmul(J31L,PDstandardNth13gt13)),kmadd(dJ111L,PDstandardNth1gt13,kmadd(ToReal(2),kmul(J21L,kmul(J31L,PDstandardNth23gt13)),kmadd(dJ211L,PDstandardNth2gt13,kmadd(dJ311L,PDstandardNth3gt13,kmadd(PDstandardNth11gt13,kmul(J11L,J11L),kmadd(PDstandardNth22gt13,kmul(J21L,J21L),kmul(PDstandardNth33gt13,kmul(J31L,J31L))))))))));
      
      JacPDstandardNth11gt22 = 
        kmadd(ToReal(2),kmul(J11L,kmul(J21L,PDstandardNth12gt22)),kmadd(ToReal(2),kmul(J11L,kmul(J31L,PDstandardNth13gt22)),kmadd(dJ111L,PDstandardNth1gt22,kmadd(ToReal(2),kmul(J21L,kmul(J31L,PDstandardNth23gt22)),kmadd(dJ211L,PDstandardNth2gt22,kmadd(dJ311L,PDstandardNth3gt22,kmadd(PDstandardNth11gt22,kmul(J11L,J11L),kmadd(PDstandardNth22gt22,kmul(J21L,J21L),kmul(PDstandardNth33gt22,kmul(J31L,J31L))))))))));
      
      JacPDstandardNth11gt23 = 
        kmadd(ToReal(2),kmul(J11L,kmul(J21L,PDstandardNth12gt23)),kmadd(ToReal(2),kmul(J11L,kmul(J31L,PDstandardNth13gt23)),kmadd(dJ111L,PDstandardNth1gt23,kmadd(ToReal(2),kmul(J21L,kmul(J31L,PDstandardNth23gt23)),kmadd(dJ211L,PDstandardNth2gt23,kmadd(dJ311L,PDstandardNth3gt23,kmadd(PDstandardNth11gt23,kmul(J11L,J11L),kmadd(PDstandardNth22gt23,kmul(J21L,J21L),kmul(PDstandardNth33gt23,kmul(J31L,J31L))))))))));
      
      JacPDstandardNth11gt33 = 
        kmadd(ToReal(2),kmul(J11L,kmul(J21L,PDstandardNth12gt33)),kmadd(ToReal(2),kmul(J11L,kmul(J31L,PDstandardNth13gt33)),kmadd(dJ111L,PDstandardNth1gt33,kmadd(ToReal(2),kmul(J21L,kmul(J31L,PDstandardNth23gt33)),kmadd(dJ211L,PDstandardNth2gt33,kmadd(dJ311L,PDstandardNth3gt33,kmadd(PDstandardNth11gt33,kmul(J11L,J11L),kmadd(PDstandardNth22gt33,kmul(J21L,J21L),kmul(PDstandardNth33gt33,kmul(J31L,J31L))))))))));
      
      JacPDstandardNth11phi = 
        kmadd(ToReal(2),kmul(J11L,kmul(J21L,PDstandardNth12phi)),kmadd(ToReal(2),kmul(J11L,kmul(J31L,PDstandardNth13phi)),kmadd(dJ111L,PDstandardNth1phi,kmadd(ToReal(2),kmul(J21L,kmul(J31L,PDstandardNth23phi)),kmadd(dJ211L,PDstandardNth2phi,kmadd(dJ311L,PDstandardNth3phi,kmadd(PDstandardNth11phi,kmul(J11L,J11L),kmadd(PDstandardNth22phi,kmul(J21L,J21L),kmul(PDstandardNth33phi,kmul(J31L,J31L))))))))));
      
      JacPDstandardNth22alpha = 
        kmadd(ToReal(2),kmul(J12L,kmul(J22L,PDstandardNth12alpha)),kmadd(ToReal(2),kmul(J12L,kmul(J32L,PDstandardNth13alpha)),kmadd(dJ122L,PDstandardNth1alpha,kmadd(ToReal(2),kmul(J22L,kmul(J32L,PDstandardNth23alpha)),kmadd(dJ222L,PDstandardNth2alpha,kmadd(dJ322L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,kmul(J12L,J12L),kmadd(PDstandardNth22alpha,kmul(J22L,J22L),kmul(PDstandardNth33alpha,kmul(J32L,J32L))))))))));
      
      JacPDstandardNth22gt11 = 
        kmadd(ToReal(2),kmul(J12L,kmul(J22L,PDstandardNth12gt11)),kmadd(ToReal(2),kmul(J12L,kmul(J32L,PDstandardNth13gt11)),kmadd(dJ122L,PDstandardNth1gt11,kmadd(ToReal(2),kmul(J22L,kmul(J32L,PDstandardNth23gt11)),kmadd(dJ222L,PDstandardNth2gt11,kmadd(dJ322L,PDstandardNth3gt11,kmadd(PDstandardNth11gt11,kmul(J12L,J12L),kmadd(PDstandardNth22gt11,kmul(J22L,J22L),kmul(PDstandardNth33gt11,kmul(J32L,J32L))))))))));
      
      JacPDstandardNth22gt12 = 
        kmadd(ToReal(2),kmul(J12L,kmul(J22L,PDstandardNth12gt12)),kmadd(ToReal(2),kmul(J12L,kmul(J32L,PDstandardNth13gt12)),kmadd(dJ122L,PDstandardNth1gt12,kmadd(ToReal(2),kmul(J22L,kmul(J32L,PDstandardNth23gt12)),kmadd(dJ222L,PDstandardNth2gt12,kmadd(dJ322L,PDstandardNth3gt12,kmadd(PDstandardNth11gt12,kmul(J12L,J12L),kmadd(PDstandardNth22gt12,kmul(J22L,J22L),kmul(PDstandardNth33gt12,kmul(J32L,J32L))))))))));
      
      JacPDstandardNth22gt13 = 
        kmadd(ToReal(2),kmul(J12L,kmul(J22L,PDstandardNth12gt13)),kmadd(ToReal(2),kmul(J12L,kmul(J32L,PDstandardNth13gt13)),kmadd(dJ122L,PDstandardNth1gt13,kmadd(ToReal(2),kmul(J22L,kmul(J32L,PDstandardNth23gt13)),kmadd(dJ222L,PDstandardNth2gt13,kmadd(dJ322L,PDstandardNth3gt13,kmadd(PDstandardNth11gt13,kmul(J12L,J12L),kmadd(PDstandardNth22gt13,kmul(J22L,J22L),kmul(PDstandardNth33gt13,kmul(J32L,J32L))))))))));
      
      JacPDstandardNth22gt22 = 
        kmadd(ToReal(2),kmul(J12L,kmul(J22L,PDstandardNth12gt22)),kmadd(ToReal(2),kmul(J12L,kmul(J32L,PDstandardNth13gt22)),kmadd(dJ122L,PDstandardNth1gt22,kmadd(ToReal(2),kmul(J22L,kmul(J32L,PDstandardNth23gt22)),kmadd(dJ222L,PDstandardNth2gt22,kmadd(dJ322L,PDstandardNth3gt22,kmadd(PDstandardNth11gt22,kmul(J12L,J12L),kmadd(PDstandardNth22gt22,kmul(J22L,J22L),kmul(PDstandardNth33gt22,kmul(J32L,J32L))))))))));
      
      JacPDstandardNth22gt23 = 
        kmadd(ToReal(2),kmul(J12L,kmul(J22L,PDstandardNth12gt23)),kmadd(ToReal(2),kmul(J12L,kmul(J32L,PDstandardNth13gt23)),kmadd(dJ122L,PDstandardNth1gt23,kmadd(ToReal(2),kmul(J22L,kmul(J32L,PDstandardNth23gt23)),kmadd(dJ222L,PDstandardNth2gt23,kmadd(dJ322L,PDstandardNth3gt23,kmadd(PDstandardNth11gt23,kmul(J12L,J12L),kmadd(PDstandardNth22gt23,kmul(J22L,J22L),kmul(PDstandardNth33gt23,kmul(J32L,J32L))))))))));
      
      JacPDstandardNth22gt33 = 
        kmadd(ToReal(2),kmul(J12L,kmul(J22L,PDstandardNth12gt33)),kmadd(ToReal(2),kmul(J12L,kmul(J32L,PDstandardNth13gt33)),kmadd(dJ122L,PDstandardNth1gt33,kmadd(ToReal(2),kmul(J22L,kmul(J32L,PDstandardNth23gt33)),kmadd(dJ222L,PDstandardNth2gt33,kmadd(dJ322L,PDstandardNth3gt33,kmadd(PDstandardNth11gt33,kmul(J12L,J12L),kmadd(PDstandardNth22gt33,kmul(J22L,J22L),kmul(PDstandardNth33gt33,kmul(J32L,J32L))))))))));
      
      JacPDstandardNth22phi = 
        kmadd(ToReal(2),kmul(J12L,kmul(J22L,PDstandardNth12phi)),kmadd(ToReal(2),kmul(J12L,kmul(J32L,PDstandardNth13phi)),kmadd(dJ122L,PDstandardNth1phi,kmadd(ToReal(2),kmul(J22L,kmul(J32L,PDstandardNth23phi)),kmadd(dJ222L,PDstandardNth2phi,kmadd(dJ322L,PDstandardNth3phi,kmadd(PDstandardNth11phi,kmul(J12L,J12L),kmadd(PDstandardNth22phi,kmul(J22L,J22L),kmul(PDstandardNth33phi,kmul(J32L,J32L))))))))));
      
      JacPDstandardNth33alpha = 
        kmadd(ToReal(2),kmul(J13L,kmul(J23L,PDstandardNth12alpha)),kmadd(ToReal(2),kmul(J13L,kmul(J33L,PDstandardNth13alpha)),kmadd(dJ133L,PDstandardNth1alpha,kmadd(ToReal(2),kmul(J23L,kmul(J33L,PDstandardNth23alpha)),kmadd(dJ233L,PDstandardNth2alpha,kmadd(dJ333L,PDstandardNth3alpha,kmadd(PDstandardNth11alpha,kmul(J13L,J13L),kmadd(PDstandardNth22alpha,kmul(J23L,J23L),kmul(PDstandardNth33alpha,kmul(J33L,J33L))))))))));
      
      JacPDstandardNth33gt11 = 
        kmadd(ToReal(2),kmul(J13L,kmul(J23L,PDstandardNth12gt11)),kmadd(ToReal(2),kmul(J13L,kmul(J33L,PDstandardNth13gt11)),kmadd(dJ133L,PDstandardNth1gt11,kmadd(ToReal(2),kmul(J23L,kmul(J33L,PDstandardNth23gt11)),kmadd(dJ233L,PDstandardNth2gt11,kmadd(dJ333L,PDstandardNth3gt11,kmadd(PDstandardNth11gt11,kmul(J13L,J13L),kmadd(PDstandardNth22gt11,kmul(J23L,J23L),kmul(PDstandardNth33gt11,kmul(J33L,J33L))))))))));
      
      JacPDstandardNth33gt12 = 
        kmadd(ToReal(2),kmul(J13L,kmul(J23L,PDstandardNth12gt12)),kmadd(ToReal(2),kmul(J13L,kmul(J33L,PDstandardNth13gt12)),kmadd(dJ133L,PDstandardNth1gt12,kmadd(ToReal(2),kmul(J23L,kmul(J33L,PDstandardNth23gt12)),kmadd(dJ233L,PDstandardNth2gt12,kmadd(dJ333L,PDstandardNth3gt12,kmadd(PDstandardNth11gt12,kmul(J13L,J13L),kmadd(PDstandardNth22gt12,kmul(J23L,J23L),kmul(PDstandardNth33gt12,kmul(J33L,J33L))))))))));
      
      JacPDstandardNth33gt13 = 
        kmadd(ToReal(2),kmul(J13L,kmul(J23L,PDstandardNth12gt13)),kmadd(ToReal(2),kmul(J13L,kmul(J33L,PDstandardNth13gt13)),kmadd(dJ133L,PDstandardNth1gt13,kmadd(ToReal(2),kmul(J23L,kmul(J33L,PDstandardNth23gt13)),kmadd(dJ233L,PDstandardNth2gt13,kmadd(dJ333L,PDstandardNth3gt13,kmadd(PDstandardNth11gt13,kmul(J13L,J13L),kmadd(PDstandardNth22gt13,kmul(J23L,J23L),kmul(PDstandardNth33gt13,kmul(J33L,J33L))))))))));
      
      JacPDstandardNth33gt22 = 
        kmadd(ToReal(2),kmul(J13L,kmul(J23L,PDstandardNth12gt22)),kmadd(ToReal(2),kmul(J13L,kmul(J33L,PDstandardNth13gt22)),kmadd(dJ133L,PDstandardNth1gt22,kmadd(ToReal(2),kmul(J23L,kmul(J33L,PDstandardNth23gt22)),kmadd(dJ233L,PDstandardNth2gt22,kmadd(dJ333L,PDstandardNth3gt22,kmadd(PDstandardNth11gt22,kmul(J13L,J13L),kmadd(PDstandardNth22gt22,kmul(J23L,J23L),kmul(PDstandardNth33gt22,kmul(J33L,J33L))))))))));
      
      JacPDstandardNth33gt23 = 
        kmadd(ToReal(2),kmul(J13L,kmul(J23L,PDstandardNth12gt23)),kmadd(ToReal(2),kmul(J13L,kmul(J33L,PDstandardNth13gt23)),kmadd(dJ133L,PDstandardNth1gt23,kmadd(ToReal(2),kmul(J23L,kmul(J33L,PDstandardNth23gt23)),kmadd(dJ233L,PDstandardNth2gt23,kmadd(dJ333L,PDstandardNth3gt23,kmadd(PDstandardNth11gt23,kmul(J13L,J13L),kmadd(PDstandardNth22gt23,kmul(J23L,J23L),kmul(PDstandardNth33gt23,kmul(J33L,J33L))))))))));
      
      JacPDstandardNth33gt33 = 
        kmadd(ToReal(2),kmul(J13L,kmul(J23L,PDstandardNth12gt33)),kmadd(ToReal(2),kmul(J13L,kmul(J33L,PDstandardNth13gt33)),kmadd(dJ133L,PDstandardNth1gt33,kmadd(ToReal(2),kmul(J23L,kmul(J33L,PDstandardNth23gt33)),kmadd(dJ233L,PDstandardNth2gt33,kmadd(dJ333L,PDstandardNth3gt33,kmadd(PDstandardNth11gt33,kmul(J13L,J13L),kmadd(PDstandardNth22gt33,kmul(J23L,J23L),kmul(PDstandardNth33gt33,kmul(J33L,J33L))))))))));
      
      JacPDstandardNth33phi = 
        kmadd(ToReal(2),kmul(J13L,kmul(J23L,PDstandardNth12phi)),kmadd(ToReal(2),kmul(J13L,kmul(J33L,PDstandardNth13phi)),kmadd(dJ133L,PDstandardNth1phi,kmadd(ToReal(2),kmul(J23L,kmul(J33L,PDstandardNth23phi)),kmadd(dJ233L,PDstandardNth2phi,kmadd(dJ333L,PDstandardNth3phi,kmadd(PDstandardNth11phi,kmul(J13L,J13L),kmadd(PDstandardNth22phi,kmul(J23L,J23L),kmul(PDstandardNth33phi,kmul(J33L,J33L))))))))));
      
      JacPDstandardNth12alpha = 
        kmadd(J11L,kmul(J12L,PDstandardNth11alpha),kmadd(J12L,kmul(J21L,PDstandardNth12alpha),kmadd(J11L,kmul(J22L,PDstandardNth12alpha),kmadd(J12L,kmul(J31L,PDstandardNth13alpha),kmadd(J11L,kmul(J32L,PDstandardNth13alpha),kmadd(dJ112L,PDstandardNth1alpha,kmadd(J21L,kmul(J22L,PDstandardNth22alpha),kmadd(J22L,kmul(J31L,PDstandardNth23alpha),kmadd(J21L,kmul(J32L,PDstandardNth23alpha),kmadd(dJ212L,PDstandardNth2alpha,kmadd(J31L,kmul(J32L,PDstandardNth33alpha),kmul(dJ312L,PDstandardNth3alpha))))))))))));
      
      JacPDstandardNth12gt11 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11gt11),kmadd(J12L,kmul(J21L,PDstandardNth12gt11),kmadd(J11L,kmul(J22L,PDstandardNth12gt11),kmadd(J12L,kmul(J31L,PDstandardNth13gt11),kmadd(J11L,kmul(J32L,PDstandardNth13gt11),kmadd(dJ112L,PDstandardNth1gt11,kmadd(J21L,kmul(J22L,PDstandardNth22gt11),kmadd(J22L,kmul(J31L,PDstandardNth23gt11),kmadd(J21L,kmul(J32L,PDstandardNth23gt11),kmadd(dJ212L,PDstandardNth2gt11,kmadd(J31L,kmul(J32L,PDstandardNth33gt11),kmul(dJ312L,PDstandardNth3gt11))))))))))));
      
      JacPDstandardNth12gt12 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11gt12),kmadd(J12L,kmul(J21L,PDstandardNth12gt12),kmadd(J11L,kmul(J22L,PDstandardNth12gt12),kmadd(J12L,kmul(J31L,PDstandardNth13gt12),kmadd(J11L,kmul(J32L,PDstandardNth13gt12),kmadd(dJ112L,PDstandardNth1gt12,kmadd(J21L,kmul(J22L,PDstandardNth22gt12),kmadd(J22L,kmul(J31L,PDstandardNth23gt12),kmadd(J21L,kmul(J32L,PDstandardNth23gt12),kmadd(dJ212L,PDstandardNth2gt12,kmadd(J31L,kmul(J32L,PDstandardNth33gt12),kmul(dJ312L,PDstandardNth3gt12))))))))))));
      
      JacPDstandardNth12gt13 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11gt13),kmadd(J12L,kmul(J21L,PDstandardNth12gt13),kmadd(J11L,kmul(J22L,PDstandardNth12gt13),kmadd(J12L,kmul(J31L,PDstandardNth13gt13),kmadd(J11L,kmul(J32L,PDstandardNth13gt13),kmadd(dJ112L,PDstandardNth1gt13,kmadd(J21L,kmul(J22L,PDstandardNth22gt13),kmadd(J22L,kmul(J31L,PDstandardNth23gt13),kmadd(J21L,kmul(J32L,PDstandardNth23gt13),kmadd(dJ212L,PDstandardNth2gt13,kmadd(J31L,kmul(J32L,PDstandardNth33gt13),kmul(dJ312L,PDstandardNth3gt13))))))))))));
      
      JacPDstandardNth12gt22 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11gt22),kmadd(J12L,kmul(J21L,PDstandardNth12gt22),kmadd(J11L,kmul(J22L,PDstandardNth12gt22),kmadd(J12L,kmul(J31L,PDstandardNth13gt22),kmadd(J11L,kmul(J32L,PDstandardNth13gt22),kmadd(dJ112L,PDstandardNth1gt22,kmadd(J21L,kmul(J22L,PDstandardNth22gt22),kmadd(J22L,kmul(J31L,PDstandardNth23gt22),kmadd(J21L,kmul(J32L,PDstandardNth23gt22),kmadd(dJ212L,PDstandardNth2gt22,kmadd(J31L,kmul(J32L,PDstandardNth33gt22),kmul(dJ312L,PDstandardNth3gt22))))))))))));
      
      JacPDstandardNth12gt23 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11gt23),kmadd(J12L,kmul(J21L,PDstandardNth12gt23),kmadd(J11L,kmul(J22L,PDstandardNth12gt23),kmadd(J12L,kmul(J31L,PDstandardNth13gt23),kmadd(J11L,kmul(J32L,PDstandardNth13gt23),kmadd(dJ112L,PDstandardNth1gt23,kmadd(J21L,kmul(J22L,PDstandardNth22gt23),kmadd(J22L,kmul(J31L,PDstandardNth23gt23),kmadd(J21L,kmul(J32L,PDstandardNth23gt23),kmadd(dJ212L,PDstandardNth2gt23,kmadd(J31L,kmul(J32L,PDstandardNth33gt23),kmul(dJ312L,PDstandardNth3gt23))))))))))));
      
      JacPDstandardNth12gt33 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11gt33),kmadd(J12L,kmul(J21L,PDstandardNth12gt33),kmadd(J11L,kmul(J22L,PDstandardNth12gt33),kmadd(J12L,kmul(J31L,PDstandardNth13gt33),kmadd(J11L,kmul(J32L,PDstandardNth13gt33),kmadd(dJ112L,PDstandardNth1gt33,kmadd(J21L,kmul(J22L,PDstandardNth22gt33),kmadd(J22L,kmul(J31L,PDstandardNth23gt33),kmadd(J21L,kmul(J32L,PDstandardNth23gt33),kmadd(dJ212L,PDstandardNth2gt33,kmadd(J31L,kmul(J32L,PDstandardNth33gt33),kmul(dJ312L,PDstandardNth3gt33))))))))))));
      
      JacPDstandardNth12phi = 
        kmadd(J11L,kmul(J12L,PDstandardNth11phi),kmadd(J12L,kmul(J21L,PDstandardNth12phi),kmadd(J11L,kmul(J22L,PDstandardNth12phi),kmadd(J12L,kmul(J31L,PDstandardNth13phi),kmadd(J11L,kmul(J32L,PDstandardNth13phi),kmadd(dJ112L,PDstandardNth1phi,kmadd(J21L,kmul(J22L,PDstandardNth22phi),kmadd(J22L,kmul(J31L,PDstandardNth23phi),kmadd(J21L,kmul(J32L,PDstandardNth23phi),kmadd(dJ212L,PDstandardNth2phi,kmadd(J31L,kmul(J32L,PDstandardNth33phi),kmul(dJ312L,PDstandardNth3phi))))))))))));
      
      JacPDstandardNth13alpha = 
        kmadd(J11L,kmul(J13L,PDstandardNth11alpha),kmadd(J13L,kmul(J21L,PDstandardNth12alpha),kmadd(J11L,kmul(J23L,PDstandardNth12alpha),kmadd(J13L,kmul(J31L,PDstandardNth13alpha),kmadd(J11L,kmul(J33L,PDstandardNth13alpha),kmadd(dJ113L,PDstandardNth1alpha,kmadd(J21L,kmul(J23L,PDstandardNth22alpha),kmadd(J23L,kmul(J31L,PDstandardNth23alpha),kmadd(J21L,kmul(J33L,PDstandardNth23alpha),kmadd(dJ213L,PDstandardNth2alpha,kmadd(J31L,kmul(J33L,PDstandardNth33alpha),kmul(dJ313L,PDstandardNth3alpha))))))))))));
      
      JacPDstandardNth13gt11 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11gt11),kmadd(J13L,kmul(J21L,PDstandardNth12gt11),kmadd(J11L,kmul(J23L,PDstandardNth12gt11),kmadd(J13L,kmul(J31L,PDstandardNth13gt11),kmadd(J11L,kmul(J33L,PDstandardNth13gt11),kmadd(dJ113L,PDstandardNth1gt11,kmadd(J21L,kmul(J23L,PDstandardNth22gt11),kmadd(J23L,kmul(J31L,PDstandardNth23gt11),kmadd(J21L,kmul(J33L,PDstandardNth23gt11),kmadd(dJ213L,PDstandardNth2gt11,kmadd(J31L,kmul(J33L,PDstandardNth33gt11),kmul(dJ313L,PDstandardNth3gt11))))))))))));
      
      JacPDstandardNth13gt12 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11gt12),kmadd(J13L,kmul(J21L,PDstandardNth12gt12),kmadd(J11L,kmul(J23L,PDstandardNth12gt12),kmadd(J13L,kmul(J31L,PDstandardNth13gt12),kmadd(J11L,kmul(J33L,PDstandardNth13gt12),kmadd(dJ113L,PDstandardNth1gt12,kmadd(J21L,kmul(J23L,PDstandardNth22gt12),kmadd(J23L,kmul(J31L,PDstandardNth23gt12),kmadd(J21L,kmul(J33L,PDstandardNth23gt12),kmadd(dJ213L,PDstandardNth2gt12,kmadd(J31L,kmul(J33L,PDstandardNth33gt12),kmul(dJ313L,PDstandardNth3gt12))))))))))));
      
      JacPDstandardNth13gt13 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11gt13),kmadd(J13L,kmul(J21L,PDstandardNth12gt13),kmadd(J11L,kmul(J23L,PDstandardNth12gt13),kmadd(J13L,kmul(J31L,PDstandardNth13gt13),kmadd(J11L,kmul(J33L,PDstandardNth13gt13),kmadd(dJ113L,PDstandardNth1gt13,kmadd(J21L,kmul(J23L,PDstandardNth22gt13),kmadd(J23L,kmul(J31L,PDstandardNth23gt13),kmadd(J21L,kmul(J33L,PDstandardNth23gt13),kmadd(dJ213L,PDstandardNth2gt13,kmadd(J31L,kmul(J33L,PDstandardNth33gt13),kmul(dJ313L,PDstandardNth3gt13))))))))))));
      
      JacPDstandardNth13gt22 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11gt22),kmadd(J13L,kmul(J21L,PDstandardNth12gt22),kmadd(J11L,kmul(J23L,PDstandardNth12gt22),kmadd(J13L,kmul(J31L,PDstandardNth13gt22),kmadd(J11L,kmul(J33L,PDstandardNth13gt22),kmadd(dJ113L,PDstandardNth1gt22,kmadd(J21L,kmul(J23L,PDstandardNth22gt22),kmadd(J23L,kmul(J31L,PDstandardNth23gt22),kmadd(J21L,kmul(J33L,PDstandardNth23gt22),kmadd(dJ213L,PDstandardNth2gt22,kmadd(J31L,kmul(J33L,PDstandardNth33gt22),kmul(dJ313L,PDstandardNth3gt22))))))))))));
      
      JacPDstandardNth13gt23 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11gt23),kmadd(J13L,kmul(J21L,PDstandardNth12gt23),kmadd(J11L,kmul(J23L,PDstandardNth12gt23),kmadd(J13L,kmul(J31L,PDstandardNth13gt23),kmadd(J11L,kmul(J33L,PDstandardNth13gt23),kmadd(dJ113L,PDstandardNth1gt23,kmadd(J21L,kmul(J23L,PDstandardNth22gt23),kmadd(J23L,kmul(J31L,PDstandardNth23gt23),kmadd(J21L,kmul(J33L,PDstandardNth23gt23),kmadd(dJ213L,PDstandardNth2gt23,kmadd(J31L,kmul(J33L,PDstandardNth33gt23),kmul(dJ313L,PDstandardNth3gt23))))))))))));
      
      JacPDstandardNth13gt33 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11gt33),kmadd(J13L,kmul(J21L,PDstandardNth12gt33),kmadd(J11L,kmul(J23L,PDstandardNth12gt33),kmadd(J13L,kmul(J31L,PDstandardNth13gt33),kmadd(J11L,kmul(J33L,PDstandardNth13gt33),kmadd(dJ113L,PDstandardNth1gt33,kmadd(J21L,kmul(J23L,PDstandardNth22gt33),kmadd(J23L,kmul(J31L,PDstandardNth23gt33),kmadd(J21L,kmul(J33L,PDstandardNth23gt33),kmadd(dJ213L,PDstandardNth2gt33,kmadd(J31L,kmul(J33L,PDstandardNth33gt33),kmul(dJ313L,PDstandardNth3gt33))))))))))));
      
      JacPDstandardNth13phi = 
        kmadd(J11L,kmul(J13L,PDstandardNth11phi),kmadd(J13L,kmul(J21L,PDstandardNth12phi),kmadd(J11L,kmul(J23L,PDstandardNth12phi),kmadd(J13L,kmul(J31L,PDstandardNth13phi),kmadd(J11L,kmul(J33L,PDstandardNth13phi),kmadd(dJ113L,PDstandardNth1phi,kmadd(J21L,kmul(J23L,PDstandardNth22phi),kmadd(J23L,kmul(J31L,PDstandardNth23phi),kmadd(J21L,kmul(J33L,PDstandardNth23phi),kmadd(dJ213L,PDstandardNth2phi,kmadd(J31L,kmul(J33L,PDstandardNth33phi),kmul(dJ313L,PDstandardNth3phi))))))))))));
      
      JacPDstandardNth21gt11 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11gt11),kmadd(J12L,kmul(J21L,PDstandardNth12gt11),kmadd(J11L,kmul(J22L,PDstandardNth12gt11),kmadd(J12L,kmul(J31L,PDstandardNth13gt11),kmadd(J11L,kmul(J32L,PDstandardNth13gt11),kmadd(dJ112L,PDstandardNth1gt11,kmadd(J21L,kmul(J22L,PDstandardNth22gt11),kmadd(J22L,kmul(J31L,PDstandardNth23gt11),kmadd(J21L,kmul(J32L,PDstandardNth23gt11),kmadd(dJ212L,PDstandardNth2gt11,kmadd(J31L,kmul(J32L,PDstandardNth33gt11),kmul(dJ312L,PDstandardNth3gt11))))))))))));
      
      JacPDstandardNth21gt12 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11gt12),kmadd(J12L,kmul(J21L,PDstandardNth12gt12),kmadd(J11L,kmul(J22L,PDstandardNth12gt12),kmadd(J12L,kmul(J31L,PDstandardNth13gt12),kmadd(J11L,kmul(J32L,PDstandardNth13gt12),kmadd(dJ112L,PDstandardNth1gt12,kmadd(J21L,kmul(J22L,PDstandardNth22gt12),kmadd(J22L,kmul(J31L,PDstandardNth23gt12),kmadd(J21L,kmul(J32L,PDstandardNth23gt12),kmadd(dJ212L,PDstandardNth2gt12,kmadd(J31L,kmul(J32L,PDstandardNth33gt12),kmul(dJ312L,PDstandardNth3gt12))))))))))));
      
      JacPDstandardNth21gt13 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11gt13),kmadd(J12L,kmul(J21L,PDstandardNth12gt13),kmadd(J11L,kmul(J22L,PDstandardNth12gt13),kmadd(J12L,kmul(J31L,PDstandardNth13gt13),kmadd(J11L,kmul(J32L,PDstandardNth13gt13),kmadd(dJ112L,PDstandardNth1gt13,kmadd(J21L,kmul(J22L,PDstandardNth22gt13),kmadd(J22L,kmul(J31L,PDstandardNth23gt13),kmadd(J21L,kmul(J32L,PDstandardNth23gt13),kmadd(dJ212L,PDstandardNth2gt13,kmadd(J31L,kmul(J32L,PDstandardNth33gt13),kmul(dJ312L,PDstandardNth3gt13))))))))))));
      
      JacPDstandardNth21gt22 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11gt22),kmadd(J12L,kmul(J21L,PDstandardNth12gt22),kmadd(J11L,kmul(J22L,PDstandardNth12gt22),kmadd(J12L,kmul(J31L,PDstandardNth13gt22),kmadd(J11L,kmul(J32L,PDstandardNth13gt22),kmadd(dJ112L,PDstandardNth1gt22,kmadd(J21L,kmul(J22L,PDstandardNth22gt22),kmadd(J22L,kmul(J31L,PDstandardNth23gt22),kmadd(J21L,kmul(J32L,PDstandardNth23gt22),kmadd(dJ212L,PDstandardNth2gt22,kmadd(J31L,kmul(J32L,PDstandardNth33gt22),kmul(dJ312L,PDstandardNth3gt22))))))))))));
      
      JacPDstandardNth21gt23 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11gt23),kmadd(J12L,kmul(J21L,PDstandardNth12gt23),kmadd(J11L,kmul(J22L,PDstandardNth12gt23),kmadd(J12L,kmul(J31L,PDstandardNth13gt23),kmadd(J11L,kmul(J32L,PDstandardNth13gt23),kmadd(dJ112L,PDstandardNth1gt23,kmadd(J21L,kmul(J22L,PDstandardNth22gt23),kmadd(J22L,kmul(J31L,PDstandardNth23gt23),kmadd(J21L,kmul(J32L,PDstandardNth23gt23),kmadd(dJ212L,PDstandardNth2gt23,kmadd(J31L,kmul(J32L,PDstandardNth33gt23),kmul(dJ312L,PDstandardNth3gt23))))))))))));
      
      JacPDstandardNth21gt33 = 
        kmadd(J11L,kmul(J12L,PDstandardNth11gt33),kmadd(J12L,kmul(J21L,PDstandardNth12gt33),kmadd(J11L,kmul(J22L,PDstandardNth12gt33),kmadd(J12L,kmul(J31L,PDstandardNth13gt33),kmadd(J11L,kmul(J32L,PDstandardNth13gt33),kmadd(dJ112L,PDstandardNth1gt33,kmadd(J21L,kmul(J22L,PDstandardNth22gt33),kmadd(J22L,kmul(J31L,PDstandardNth23gt33),kmadd(J21L,kmul(J32L,PDstandardNth23gt33),kmadd(dJ212L,PDstandardNth2gt33,kmadd(J31L,kmul(J32L,PDstandardNth33gt33),kmul(dJ312L,PDstandardNth3gt33))))))))))));
      
      JacPDstandardNth23alpha = 
        kmadd(J12L,kmul(J13L,PDstandardNth11alpha),kmadd(J13L,kmul(J22L,PDstandardNth12alpha),kmadd(J12L,kmul(J23L,PDstandardNth12alpha),kmadd(J13L,kmul(J32L,PDstandardNth13alpha),kmadd(J12L,kmul(J33L,PDstandardNth13alpha),kmadd(dJ123L,PDstandardNth1alpha,kmadd(J22L,kmul(J23L,PDstandardNth22alpha),kmadd(J23L,kmul(J32L,PDstandardNth23alpha),kmadd(J22L,kmul(J33L,PDstandardNth23alpha),kmadd(dJ223L,PDstandardNth2alpha,kmadd(J32L,kmul(J33L,PDstandardNth33alpha),kmul(dJ323L,PDstandardNth3alpha))))))))))));
      
      JacPDstandardNth23gt11 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11gt11),kmadd(J13L,kmul(J22L,PDstandardNth12gt11),kmadd(J12L,kmul(J23L,PDstandardNth12gt11),kmadd(J13L,kmul(J32L,PDstandardNth13gt11),kmadd(J12L,kmul(J33L,PDstandardNth13gt11),kmadd(dJ123L,PDstandardNth1gt11,kmadd(J22L,kmul(J23L,PDstandardNth22gt11),kmadd(J23L,kmul(J32L,PDstandardNth23gt11),kmadd(J22L,kmul(J33L,PDstandardNth23gt11),kmadd(dJ223L,PDstandardNth2gt11,kmadd(J32L,kmul(J33L,PDstandardNth33gt11),kmul(dJ323L,PDstandardNth3gt11))))))))))));
      
      JacPDstandardNth23gt12 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11gt12),kmadd(J13L,kmul(J22L,PDstandardNth12gt12),kmadd(J12L,kmul(J23L,PDstandardNth12gt12),kmadd(J13L,kmul(J32L,PDstandardNth13gt12),kmadd(J12L,kmul(J33L,PDstandardNth13gt12),kmadd(dJ123L,PDstandardNth1gt12,kmadd(J22L,kmul(J23L,PDstandardNth22gt12),kmadd(J23L,kmul(J32L,PDstandardNth23gt12),kmadd(J22L,kmul(J33L,PDstandardNth23gt12),kmadd(dJ223L,PDstandardNth2gt12,kmadd(J32L,kmul(J33L,PDstandardNth33gt12),kmul(dJ323L,PDstandardNth3gt12))))))))))));
      
      JacPDstandardNth23gt13 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11gt13),kmadd(J13L,kmul(J22L,PDstandardNth12gt13),kmadd(J12L,kmul(J23L,PDstandardNth12gt13),kmadd(J13L,kmul(J32L,PDstandardNth13gt13),kmadd(J12L,kmul(J33L,PDstandardNth13gt13),kmadd(dJ123L,PDstandardNth1gt13,kmadd(J22L,kmul(J23L,PDstandardNth22gt13),kmadd(J23L,kmul(J32L,PDstandardNth23gt13),kmadd(J22L,kmul(J33L,PDstandardNth23gt13),kmadd(dJ223L,PDstandardNth2gt13,kmadd(J32L,kmul(J33L,PDstandardNth33gt13),kmul(dJ323L,PDstandardNth3gt13))))))))))));
      
      JacPDstandardNth23gt22 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11gt22),kmadd(J13L,kmul(J22L,PDstandardNth12gt22),kmadd(J12L,kmul(J23L,PDstandardNth12gt22),kmadd(J13L,kmul(J32L,PDstandardNth13gt22),kmadd(J12L,kmul(J33L,PDstandardNth13gt22),kmadd(dJ123L,PDstandardNth1gt22,kmadd(J22L,kmul(J23L,PDstandardNth22gt22),kmadd(J23L,kmul(J32L,PDstandardNth23gt22),kmadd(J22L,kmul(J33L,PDstandardNth23gt22),kmadd(dJ223L,PDstandardNth2gt22,kmadd(J32L,kmul(J33L,PDstandardNth33gt22),kmul(dJ323L,PDstandardNth3gt22))))))))))));
      
      JacPDstandardNth23gt23 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11gt23),kmadd(J13L,kmul(J22L,PDstandardNth12gt23),kmadd(J12L,kmul(J23L,PDstandardNth12gt23),kmadd(J13L,kmul(J32L,PDstandardNth13gt23),kmadd(J12L,kmul(J33L,PDstandardNth13gt23),kmadd(dJ123L,PDstandardNth1gt23,kmadd(J22L,kmul(J23L,PDstandardNth22gt23),kmadd(J23L,kmul(J32L,PDstandardNth23gt23),kmadd(J22L,kmul(J33L,PDstandardNth23gt23),kmadd(dJ223L,PDstandardNth2gt23,kmadd(J32L,kmul(J33L,PDstandardNth33gt23),kmul(dJ323L,PDstandardNth3gt23))))))))))));
      
      JacPDstandardNth23gt33 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11gt33),kmadd(J13L,kmul(J22L,PDstandardNth12gt33),kmadd(J12L,kmul(J23L,PDstandardNth12gt33),kmadd(J13L,kmul(J32L,PDstandardNth13gt33),kmadd(J12L,kmul(J33L,PDstandardNth13gt33),kmadd(dJ123L,PDstandardNth1gt33,kmadd(J22L,kmul(J23L,PDstandardNth22gt33),kmadd(J23L,kmul(J32L,PDstandardNth23gt33),kmadd(J22L,kmul(J33L,PDstandardNth23gt33),kmadd(dJ223L,PDstandardNth2gt33,kmadd(J32L,kmul(J33L,PDstandardNth33gt33),kmul(dJ323L,PDstandardNth3gt33))))))))))));
      
      JacPDstandardNth23phi = 
        kmadd(J12L,kmul(J13L,PDstandardNth11phi),kmadd(J13L,kmul(J22L,PDstandardNth12phi),kmadd(J12L,kmul(J23L,PDstandardNth12phi),kmadd(J13L,kmul(J32L,PDstandardNth13phi),kmadd(J12L,kmul(J33L,PDstandardNth13phi),kmadd(dJ123L,PDstandardNth1phi,kmadd(J22L,kmul(J23L,PDstandardNth22phi),kmadd(J23L,kmul(J32L,PDstandardNth23phi),kmadd(J22L,kmul(J33L,PDstandardNth23phi),kmadd(dJ223L,PDstandardNth2phi,kmadd(J32L,kmul(J33L,PDstandardNth33phi),kmul(dJ323L,PDstandardNth3phi))))))))))));
      
      JacPDstandardNth31gt11 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11gt11),kmadd(J13L,kmul(J21L,PDstandardNth12gt11),kmadd(J11L,kmul(J23L,PDstandardNth12gt11),kmadd(J13L,kmul(J31L,PDstandardNth13gt11),kmadd(J11L,kmul(J33L,PDstandardNth13gt11),kmadd(dJ113L,PDstandardNth1gt11,kmadd(J21L,kmul(J23L,PDstandardNth22gt11),kmadd(J23L,kmul(J31L,PDstandardNth23gt11),kmadd(J21L,kmul(J33L,PDstandardNth23gt11),kmadd(dJ213L,PDstandardNth2gt11,kmadd(J31L,kmul(J33L,PDstandardNth33gt11),kmul(dJ313L,PDstandardNth3gt11))))))))))));
      
      JacPDstandardNth31gt12 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11gt12),kmadd(J13L,kmul(J21L,PDstandardNth12gt12),kmadd(J11L,kmul(J23L,PDstandardNth12gt12),kmadd(J13L,kmul(J31L,PDstandardNth13gt12),kmadd(J11L,kmul(J33L,PDstandardNth13gt12),kmadd(dJ113L,PDstandardNth1gt12,kmadd(J21L,kmul(J23L,PDstandardNth22gt12),kmadd(J23L,kmul(J31L,PDstandardNth23gt12),kmadd(J21L,kmul(J33L,PDstandardNth23gt12),kmadd(dJ213L,PDstandardNth2gt12,kmadd(J31L,kmul(J33L,PDstandardNth33gt12),kmul(dJ313L,PDstandardNth3gt12))))))))))));
      
      JacPDstandardNth31gt13 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11gt13),kmadd(J13L,kmul(J21L,PDstandardNth12gt13),kmadd(J11L,kmul(J23L,PDstandardNth12gt13),kmadd(J13L,kmul(J31L,PDstandardNth13gt13),kmadd(J11L,kmul(J33L,PDstandardNth13gt13),kmadd(dJ113L,PDstandardNth1gt13,kmadd(J21L,kmul(J23L,PDstandardNth22gt13),kmadd(J23L,kmul(J31L,PDstandardNth23gt13),kmadd(J21L,kmul(J33L,PDstandardNth23gt13),kmadd(dJ213L,PDstandardNth2gt13,kmadd(J31L,kmul(J33L,PDstandardNth33gt13),kmul(dJ313L,PDstandardNth3gt13))))))))))));
      
      JacPDstandardNth31gt22 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11gt22),kmadd(J13L,kmul(J21L,PDstandardNth12gt22),kmadd(J11L,kmul(J23L,PDstandardNth12gt22),kmadd(J13L,kmul(J31L,PDstandardNth13gt22),kmadd(J11L,kmul(J33L,PDstandardNth13gt22),kmadd(dJ113L,PDstandardNth1gt22,kmadd(J21L,kmul(J23L,PDstandardNth22gt22),kmadd(J23L,kmul(J31L,PDstandardNth23gt22),kmadd(J21L,kmul(J33L,PDstandardNth23gt22),kmadd(dJ213L,PDstandardNth2gt22,kmadd(J31L,kmul(J33L,PDstandardNth33gt22),kmul(dJ313L,PDstandardNth3gt22))))))))))));
      
      JacPDstandardNth31gt23 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11gt23),kmadd(J13L,kmul(J21L,PDstandardNth12gt23),kmadd(J11L,kmul(J23L,PDstandardNth12gt23),kmadd(J13L,kmul(J31L,PDstandardNth13gt23),kmadd(J11L,kmul(J33L,PDstandardNth13gt23),kmadd(dJ113L,PDstandardNth1gt23,kmadd(J21L,kmul(J23L,PDstandardNth22gt23),kmadd(J23L,kmul(J31L,PDstandardNth23gt23),kmadd(J21L,kmul(J33L,PDstandardNth23gt23),kmadd(dJ213L,PDstandardNth2gt23,kmadd(J31L,kmul(J33L,PDstandardNth33gt23),kmul(dJ313L,PDstandardNth3gt23))))))))))));
      
      JacPDstandardNth31gt33 = 
        kmadd(J11L,kmul(J13L,PDstandardNth11gt33),kmadd(J13L,kmul(J21L,PDstandardNth12gt33),kmadd(J11L,kmul(J23L,PDstandardNth12gt33),kmadd(J13L,kmul(J31L,PDstandardNth13gt33),kmadd(J11L,kmul(J33L,PDstandardNth13gt33),kmadd(dJ113L,PDstandardNth1gt33,kmadd(J21L,kmul(J23L,PDstandardNth22gt33),kmadd(J23L,kmul(J31L,PDstandardNth23gt33),kmadd(J21L,kmul(J33L,PDstandardNth23gt33),kmadd(dJ213L,PDstandardNth2gt33,kmadd(J31L,kmul(J33L,PDstandardNth33gt33),kmul(dJ313L,PDstandardNth3gt33))))))))))));
      
      JacPDstandardNth32gt11 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11gt11),kmadd(J13L,kmul(J22L,PDstandardNth12gt11),kmadd(J12L,kmul(J23L,PDstandardNth12gt11),kmadd(J13L,kmul(J32L,PDstandardNth13gt11),kmadd(J12L,kmul(J33L,PDstandardNth13gt11),kmadd(dJ123L,PDstandardNth1gt11,kmadd(J22L,kmul(J23L,PDstandardNth22gt11),kmadd(J23L,kmul(J32L,PDstandardNth23gt11),kmadd(J22L,kmul(J33L,PDstandardNth23gt11),kmadd(dJ223L,PDstandardNth2gt11,kmadd(J32L,kmul(J33L,PDstandardNth33gt11),kmul(dJ323L,PDstandardNth3gt11))))))))))));
      
      JacPDstandardNth32gt12 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11gt12),kmadd(J13L,kmul(J22L,PDstandardNth12gt12),kmadd(J12L,kmul(J23L,PDstandardNth12gt12),kmadd(J13L,kmul(J32L,PDstandardNth13gt12),kmadd(J12L,kmul(J33L,PDstandardNth13gt12),kmadd(dJ123L,PDstandardNth1gt12,kmadd(J22L,kmul(J23L,PDstandardNth22gt12),kmadd(J23L,kmul(J32L,PDstandardNth23gt12),kmadd(J22L,kmul(J33L,PDstandardNth23gt12),kmadd(dJ223L,PDstandardNth2gt12,kmadd(J32L,kmul(J33L,PDstandardNth33gt12),kmul(dJ323L,PDstandardNth3gt12))))))))))));
      
      JacPDstandardNth32gt13 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11gt13),kmadd(J13L,kmul(J22L,PDstandardNth12gt13),kmadd(J12L,kmul(J23L,PDstandardNth12gt13),kmadd(J13L,kmul(J32L,PDstandardNth13gt13),kmadd(J12L,kmul(J33L,PDstandardNth13gt13),kmadd(dJ123L,PDstandardNth1gt13,kmadd(J22L,kmul(J23L,PDstandardNth22gt13),kmadd(J23L,kmul(J32L,PDstandardNth23gt13),kmadd(J22L,kmul(J33L,PDstandardNth23gt13),kmadd(dJ223L,PDstandardNth2gt13,kmadd(J32L,kmul(J33L,PDstandardNth33gt13),kmul(dJ323L,PDstandardNth3gt13))))))))))));
      
      JacPDstandardNth32gt22 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11gt22),kmadd(J13L,kmul(J22L,PDstandardNth12gt22),kmadd(J12L,kmul(J23L,PDstandardNth12gt22),kmadd(J13L,kmul(J32L,PDstandardNth13gt22),kmadd(J12L,kmul(J33L,PDstandardNth13gt22),kmadd(dJ123L,PDstandardNth1gt22,kmadd(J22L,kmul(J23L,PDstandardNth22gt22),kmadd(J23L,kmul(J32L,PDstandardNth23gt22),kmadd(J22L,kmul(J33L,PDstandardNth23gt22),kmadd(dJ223L,PDstandardNth2gt22,kmadd(J32L,kmul(J33L,PDstandardNth33gt22),kmul(dJ323L,PDstandardNth3gt22))))))))))));
      
      JacPDstandardNth32gt23 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11gt23),kmadd(J13L,kmul(J22L,PDstandardNth12gt23),kmadd(J12L,kmul(J23L,PDstandardNth12gt23),kmadd(J13L,kmul(J32L,PDstandardNth13gt23),kmadd(J12L,kmul(J33L,PDstandardNth13gt23),kmadd(dJ123L,PDstandardNth1gt23,kmadd(J22L,kmul(J23L,PDstandardNth22gt23),kmadd(J23L,kmul(J32L,PDstandardNth23gt23),kmadd(J22L,kmul(J33L,PDstandardNth23gt23),kmadd(dJ223L,PDstandardNth2gt23,kmadd(J32L,kmul(J33L,PDstandardNth33gt23),kmul(dJ323L,PDstandardNth3gt23))))))))))));
      
      JacPDstandardNth32gt33 = 
        kmadd(J12L,kmul(J13L,PDstandardNth11gt33),kmadd(J13L,kmul(J22L,PDstandardNth12gt33),kmadd(J12L,kmul(J23L,PDstandardNth12gt33),kmadd(J13L,kmul(J32L,PDstandardNth13gt33),kmadd(J12L,kmul(J33L,PDstandardNth13gt33),kmadd(dJ123L,PDstandardNth1gt33,kmadd(J22L,kmul(J23L,PDstandardNth22gt33),kmadd(J23L,kmul(J32L,PDstandardNth23gt33),kmadd(J22L,kmul(J33L,PDstandardNth23gt33),kmadd(dJ223L,PDstandardNth2gt33,kmadd(J32L,kmul(J33L,PDstandardNth33gt33),kmul(dJ323L,PDstandardNth3gt33))))))))))));
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
      
      JacPDstandardNth1Xt1 = PDstandardNth1Xt1;
      
      JacPDstandardNth1Xt2 = PDstandardNth1Xt2;
      
      JacPDstandardNth1Xt3 = PDstandardNth1Xt3;
      
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
      
      JacPDstandardNth2Xt1 = PDstandardNth2Xt1;
      
      JacPDstandardNth2Xt2 = PDstandardNth2Xt2;
      
      JacPDstandardNth2Xt3 = PDstandardNth2Xt3;
      
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
      
      JacPDstandardNth3Xt1 = PDstandardNth3Xt1;
      
      JacPDstandardNth3Xt2 = PDstandardNth3Xt2;
      
      JacPDstandardNth3Xt3 = PDstandardNth3Xt3;
      
      JacPDstandardNth11alpha = PDstandardNth11alpha;
      
      JacPDstandardNth11gt11 = PDstandardNth11gt11;
      
      JacPDstandardNth11gt12 = PDstandardNth11gt12;
      
      JacPDstandardNth11gt13 = PDstandardNth11gt13;
      
      JacPDstandardNth11gt22 = PDstandardNth11gt22;
      
      JacPDstandardNth11gt23 = PDstandardNth11gt23;
      
      JacPDstandardNth11gt33 = PDstandardNth11gt33;
      
      JacPDstandardNth11phi = PDstandardNth11phi;
      
      JacPDstandardNth22alpha = PDstandardNth22alpha;
      
      JacPDstandardNth22gt11 = PDstandardNth22gt11;
      
      JacPDstandardNth22gt12 = PDstandardNth22gt12;
      
      JacPDstandardNth22gt13 = PDstandardNth22gt13;
      
      JacPDstandardNth22gt22 = PDstandardNth22gt22;
      
      JacPDstandardNth22gt23 = PDstandardNth22gt23;
      
      JacPDstandardNth22gt33 = PDstandardNth22gt33;
      
      JacPDstandardNth22phi = PDstandardNth22phi;
      
      JacPDstandardNth33alpha = PDstandardNth33alpha;
      
      JacPDstandardNth33gt11 = PDstandardNth33gt11;
      
      JacPDstandardNth33gt12 = PDstandardNth33gt12;
      
      JacPDstandardNth33gt13 = PDstandardNth33gt13;
      
      JacPDstandardNth33gt22 = PDstandardNth33gt22;
      
      JacPDstandardNth33gt23 = PDstandardNth33gt23;
      
      JacPDstandardNth33gt33 = PDstandardNth33gt33;
      
      JacPDstandardNth33phi = PDstandardNth33phi;
      
      JacPDstandardNth12alpha = PDstandardNth12alpha;
      
      JacPDstandardNth12gt11 = PDstandardNth12gt11;
      
      JacPDstandardNth12gt12 = PDstandardNth12gt12;
      
      JacPDstandardNth12gt13 = PDstandardNth12gt13;
      
      JacPDstandardNth12gt22 = PDstandardNth12gt22;
      
      JacPDstandardNth12gt23 = PDstandardNth12gt23;
      
      JacPDstandardNth12gt33 = PDstandardNth12gt33;
      
      JacPDstandardNth12phi = PDstandardNth12phi;
      
      JacPDstandardNth13alpha = PDstandardNth13alpha;
      
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
      
      JacPDstandardNth23alpha = PDstandardNth23alpha;
      
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
      kmul(kmsub(ToReal(2),JacPDstandardNth2gt12,JacPDstandardNth1gt22),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl123 CCTK_ATTRIBUTE_UNUSED = 
      kmul(ksub(kadd(JacPDstandardNth2gt13,JacPDstandardNth3gt12),JacPDstandardNth1gt23),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl133 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmsub(ToReal(2),JacPDstandardNth3gt13,JacPDstandardNth1gt33),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl211 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmsub(ToReal(2),JacPDstandardNth1gt12,JacPDstandardNth2gt11),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl212 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth1gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl213 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth3gt12,JacPDstandardNth2gt13)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl222 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth2gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl223 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth3gt22,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl233 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmsub(ToReal(2),JacPDstandardNth3gt23,JacPDstandardNth2gt33),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl311 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmsub(ToReal(2),JacPDstandardNth1gt13,JacPDstandardNth3gt11),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl312 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kadd(JacPDstandardNth1gt23,ksub(JacPDstandardNth2gt13,JacPDstandardNth3gt12)),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl313 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth1gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl322 CCTK_ATTRIBUTE_UNUSED = 
      kmul(kmsub(ToReal(2),JacPDstandardNth2gt23,JacPDstandardNth3gt22),ToReal(0.5));
    
    CCTK_REAL_VEC Gtl323 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth2gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtl333 CCTK_ATTRIBUTE_UNUSED = 
      kmul(JacPDstandardNth3gt33,ToReal(0.5));
    
    CCTK_REAL_VEC Gtlu111 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu11,kmadd(Gtl112,gtu12,kmul(Gtl113,gtu13)));
    
    CCTK_REAL_VEC Gtlu112 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu12,kmadd(Gtl112,gtu22,kmul(Gtl113,gtu23)));
    
    CCTK_REAL_VEC Gtlu113 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl111,gtu13,kmadd(Gtl112,gtu23,kmul(Gtl113,gtu33)));
    
    CCTK_REAL_VEC Gtlu121 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu11,kmadd(Gtl122,gtu12,kmul(Gtl123,gtu13)));
    
    CCTK_REAL_VEC Gtlu122 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu12,kmadd(Gtl122,gtu22,kmul(Gtl123,gtu23)));
    
    CCTK_REAL_VEC Gtlu123 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl112,gtu13,kmadd(Gtl122,gtu23,kmul(Gtl123,gtu33)));
    
    CCTK_REAL_VEC Gtlu131 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu11,kmadd(Gtl123,gtu12,kmul(Gtl133,gtu13)));
    
    CCTK_REAL_VEC Gtlu132 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu12,kmadd(Gtl123,gtu22,kmul(Gtl133,gtu23)));
    
    CCTK_REAL_VEC Gtlu133 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl113,gtu13,kmadd(Gtl123,gtu23,kmul(Gtl133,gtu33)));
    
    CCTK_REAL_VEC Gtlu211 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl211,gtu11,kmadd(Gtl212,gtu12,kmul(Gtl213,gtu13)));
    
    CCTK_REAL_VEC Gtlu212 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl211,gtu12,kmadd(Gtl212,gtu22,kmul(Gtl213,gtu23)));
    
    CCTK_REAL_VEC Gtlu213 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl211,gtu13,kmadd(Gtl212,gtu23,kmul(Gtl213,gtu33)));
    
    CCTK_REAL_VEC Gtlu221 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl212,gtu11,kmadd(Gtl222,gtu12,kmul(Gtl223,gtu13)));
    
    CCTK_REAL_VEC Gtlu222 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl212,gtu12,kmadd(Gtl222,gtu22,kmul(Gtl223,gtu23)));
    
    CCTK_REAL_VEC Gtlu223 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl212,gtu13,kmadd(Gtl222,gtu23,kmul(Gtl223,gtu33)));
    
    CCTK_REAL_VEC Gtlu231 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl213,gtu11,kmadd(Gtl223,gtu12,kmul(Gtl233,gtu13)));
    
    CCTK_REAL_VEC Gtlu232 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl213,gtu12,kmadd(Gtl223,gtu22,kmul(Gtl233,gtu23)));
    
    CCTK_REAL_VEC Gtlu233 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl213,gtu13,kmadd(Gtl223,gtu23,kmul(Gtl233,gtu33)));
    
    CCTK_REAL_VEC Gtlu311 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl311,gtu11,kmadd(Gtl312,gtu12,kmul(Gtl313,gtu13)));
    
    CCTK_REAL_VEC Gtlu312 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl311,gtu12,kmadd(Gtl312,gtu22,kmul(Gtl313,gtu23)));
    
    CCTK_REAL_VEC Gtlu313 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl311,gtu13,kmadd(Gtl312,gtu23,kmul(Gtl313,gtu33)));
    
    CCTK_REAL_VEC Gtlu321 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl312,gtu11,kmadd(Gtl322,gtu12,kmul(Gtl323,gtu13)));
    
    CCTK_REAL_VEC Gtlu322 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl312,gtu12,kmadd(Gtl322,gtu22,kmul(Gtl323,gtu23)));
    
    CCTK_REAL_VEC Gtlu323 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl312,gtu13,kmadd(Gtl322,gtu23,kmul(Gtl323,gtu33)));
    
    CCTK_REAL_VEC Gtlu331 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl313,gtu11,kmadd(Gtl323,gtu12,kmul(Gtl333,gtu13)));
    
    CCTK_REAL_VEC Gtlu332 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl313,gtu12,kmadd(Gtl323,gtu22,kmul(Gtl333,gtu23)));
    
    CCTK_REAL_VEC Gtlu333 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gtl313,gtu13,kmadd(Gtl323,gtu23,kmul(Gtl333,gtu33)));
    
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
      kmadd(Gt111,gtu11,kmadd(ToReal(2),kmul(Gt112,gtu12),kmadd(ToReal(2),kmul(Gt113,gtu13),kmadd(Gt122,gtu22,kmadd(ToReal(2),kmul(Gt123,gtu23),kmul(Gt133,gtu33))))));
    
    CCTK_REAL_VEC Xtn2 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt211,gtu11,kmadd(ToReal(2),kmul(Gt212,gtu12),kmadd(ToReal(2),kmul(Gt213,gtu13),kmadd(Gt222,gtu22,kmadd(ToReal(2),kmul(Gt223,gtu23),kmul(Gt233,gtu33))))));
    
    CCTK_REAL_VEC Xtn3 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt311,gtu11,kmadd(ToReal(2),kmul(Gt312,gtu12),kmadd(ToReal(2),kmul(Gt313,gtu13),kmadd(Gt322,gtu22,kmadd(ToReal(2),kmul(Gt323,gtu23),kmul(Gt333,gtu33))))));
    
    CCTK_REAL_VEC e4phi CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod == 
      1,kdiv(ToReal(1),kmul(phiL,phiL)),kexp(kmul(phiL,ToReal(4))));
    
    CCTK_REAL_VEC em4phi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),e4phi);
    
    CCTK_REAL_VEC g11 CCTK_ATTRIBUTE_UNUSED = kmul(gt11L,e4phi);
    
    CCTK_REAL_VEC g12 CCTK_ATTRIBUTE_UNUSED = kmul(gt12L,e4phi);
    
    CCTK_REAL_VEC g13 CCTK_ATTRIBUTE_UNUSED = kmul(gt13L,e4phi);
    
    CCTK_REAL_VEC g22 CCTK_ATTRIBUTE_UNUSED = kmul(gt22L,e4phi);
    
    CCTK_REAL_VEC g23 CCTK_ATTRIBUTE_UNUSED = kmul(gt23L,e4phi);
    
    CCTK_REAL_VEC g33 CCTK_ATTRIBUTE_UNUSED = kmul(gt33L,e4phi);
    
    CCTK_REAL_VEC gu11 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu11);
    
    CCTK_REAL_VEC gu12 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu12);
    
    CCTK_REAL_VEC gu13 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu13);
    
    CCTK_REAL_VEC gu22 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu22);
    
    CCTK_REAL_VEC gu23 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu23);
    
    CCTK_REAL_VEC gu33 CCTK_ATTRIBUTE_UNUSED = kmul(em4phi,gtu33);
    
    CCTK_REAL_VEC Rt11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(3),kmul(Gt111,Gtlu111),kmadd(ToReal(3),kmul(Gt112,Gtlu112),kmadd(ToReal(3),kmul(Gt113,Gtlu113),kmadd(ToReal(2),kmul(Gt211,Gtlu121),kmadd(ToReal(2),kmul(Gt212,Gtlu122),kmadd(ToReal(2),kmul(Gt213,Gtlu123),kmadd(ToReal(2),kmul(Gt311,Gtlu131),kmadd(ToReal(2),kmul(Gt312,Gtlu132),kmadd(ToReal(2),kmul(Gt313,Gtlu133),kmadd(Gt211,Gtlu211,kmadd(Gt212,Gtlu212,kmadd(Gt213,Gtlu213,kmadd(Gt311,Gtlu311,kmadd(Gt312,Gtlu312,kmadd(Gt313,Gtlu313,kmadd(gt11L,JacPDstandardNth1Xt1,kmadd(gt12L,JacPDstandardNth1Xt2,kmadd(gt13L,JacPDstandardNth1Xt3,knmsub(ToReal(0.5),kmadd(gtu11,JacPDstandardNth11gt11,kmadd(gtu12,JacPDstandardNth12gt11,kmadd(gtu13,JacPDstandardNth13gt11,kmadd(gtu12,JacPDstandardNth21gt11,kmadd(gtu22,JacPDstandardNth22gt11,kmadd(gtu23,JacPDstandardNth23gt11,kmadd(gtu13,JacPDstandardNth31gt11,kmadd(gtu33,JacPDstandardNth33gt11,kmul(gtu23,JacPDstandardNth32gt11))))))))),kmadd(Gtl111,Xtn1,kmadd(Gtl112,Xtn2,kmul(Gtl113,Xtn3))))))))))))))))))))));
    
    CCTK_REAL_VEC Rt12 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt112,Gtlu111,kmadd(Gt122,Gtlu112,kmadd(Gt123,Gtlu113,kmadd(Gt111,Gtlu121,kmadd(Gt212,Gtlu121,kmadd(Gt112,Gtlu122,kmadd(Gt222,Gtlu122,kmadd(Gt113,Gtlu123,kmadd(Gt223,Gtlu123,kmadd(Gt312,Gtlu131,kmadd(Gt322,Gtlu132,kmadd(Gt323,Gtlu133,kmadd(Gt111,Gtlu211,kmadd(Gt112,Gtlu212,kmadd(Gt113,Gtlu213,kmadd(ToReal(2),kmul(Gt211,Gtlu221),kmadd(ToReal(2),kmul(Gt212,Gtlu222),kmadd(ToReal(2),kmul(Gt213,Gtlu223),kmadd(Gt311,Gtlu231,kmadd(Gt312,Gtlu232,kmadd(Gt313,Gtlu233,kmadd(Gt311,Gtlu321,kmadd(Gt312,Gtlu322,kmadd(Gt313,Gtlu323,kmadd(ToReal(0.5),kmadd(gt12L,JacPDstandardNth1Xt1,kmadd(gt22L,JacPDstandardNth1Xt2,kmul(gt23L,JacPDstandardNth1Xt3))),kmadd(ToReal(0.5),kmadd(gt11L,JacPDstandardNth2Xt1,kmadd(gt12L,JacPDstandardNth2Xt2,kmul(gt13L,JacPDstandardNth2Xt3))),knmsub(ToReal(0.5),kmadd(gtu11,JacPDstandardNth11gt12,kmadd(gtu12,JacPDstandardNth12gt12,kmadd(gtu13,JacPDstandardNth13gt12,kmadd(gtu12,JacPDstandardNth21gt12,kmadd(gtu22,JacPDstandardNth22gt12,kmadd(gtu23,JacPDstandardNth23gt12,kmadd(gtu13,JacPDstandardNth31gt12,kmadd(gtu33,JacPDstandardNth33gt12,kmul(gtu23,JacPDstandardNth32gt12))))))))),kmadd(ToReal(0.5),kmadd(Gtl112,Xtn1,kmadd(Gtl122,Xtn2,kmul(Gtl123,Xtn3))),kmul(kmadd(Gtl211,Xtn1,kmadd(Gtl212,Xtn2,kmul(Gtl213,Xtn3))),ToReal(0.5))))))))))))))))))))))))))))));
    
    CCTK_REAL_VEC Rt13 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt113,Gtlu111,kmadd(Gt123,Gtlu112,kmadd(Gt133,Gtlu113,kmadd(Gt213,Gtlu121,kmadd(Gt223,Gtlu122,kmadd(Gt233,Gtlu123,kmadd(Gt111,Gtlu131,kmadd(Gt313,Gtlu131,kmadd(Gt112,Gtlu132,kmadd(Gt323,Gtlu132,kmadd(Gt113,Gtlu133,kmadd(Gt333,Gtlu133,kmadd(Gt211,Gtlu231,kmadd(Gt212,Gtlu232,kmadd(Gt213,Gtlu233,kmadd(Gt111,Gtlu311,kmadd(Gt112,Gtlu312,kmadd(Gt113,Gtlu313,kmadd(Gt211,Gtlu321,kmadd(Gt212,Gtlu322,kmadd(Gt213,Gtlu323,kmadd(ToReal(2),kmul(Gt311,Gtlu331),kmadd(ToReal(2),kmul(Gt312,Gtlu332),kmadd(ToReal(2),kmul(Gt313,Gtlu333),kmadd(ToReal(0.5),kmadd(gt13L,JacPDstandardNth1Xt1,kmadd(gt23L,JacPDstandardNth1Xt2,kmul(gt33L,JacPDstandardNth1Xt3))),knmsub(ToReal(0.5),kmadd(gtu11,JacPDstandardNth11gt13,kmadd(gtu12,JacPDstandardNth12gt13,kmadd(gtu13,JacPDstandardNth13gt13,kmadd(gtu12,JacPDstandardNth21gt13,kmadd(gtu22,JacPDstandardNth22gt13,kmadd(gtu23,JacPDstandardNth23gt13,kmadd(gtu13,JacPDstandardNth31gt13,kmadd(gtu33,JacPDstandardNth33gt13,kmul(gtu23,JacPDstandardNth32gt13))))))))),kmadd(ToReal(0.5),kmadd(gt11L,JacPDstandardNth3Xt1,kmadd(gt12L,JacPDstandardNth3Xt2,kmul(gt13L,JacPDstandardNth3Xt3))),kmadd(ToReal(0.5),kmadd(Gtl113,Xtn1,kmadd(Gtl123,Xtn2,kmul(Gtl133,Xtn3))),kmul(kmadd(Gtl311,Xtn1,kmadd(Gtl312,Xtn2,kmul(Gtl313,Xtn3))),ToReal(0.5))))))))))))))))))))))))))))));
    
    CCTK_REAL_VEC Rt22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt112,Gtlu121,kmadd(Gt122,Gtlu122,kmadd(Gt123,Gtlu123,kmadd(ToReal(2),kmul(Gt112,Gtlu211),kmadd(ToReal(2),kmul(Gt122,Gtlu212),kmadd(ToReal(2),kmul(Gt123,Gtlu213),kmadd(ToReal(3),kmul(Gt212,Gtlu221),kmadd(ToReal(3),kmul(Gt222,Gtlu222),kmadd(ToReal(3),kmul(Gt223,Gtlu223),kmadd(ToReal(2),kmul(Gt312,Gtlu231),kmadd(ToReal(2),kmul(Gt322,Gtlu232),kmadd(ToReal(2),kmul(Gt323,Gtlu233),kmadd(Gt312,Gtlu321,kmadd(Gt322,Gtlu322,kmadd(Gt323,Gtlu323,kmadd(gt12L,JacPDstandardNth2Xt1,kmadd(gt22L,JacPDstandardNth2Xt2,kmadd(gt23L,JacPDstandardNth2Xt3,knmsub(ToReal(0.5),kmadd(gtu11,JacPDstandardNth11gt22,kmadd(gtu12,JacPDstandardNth12gt22,kmadd(gtu13,JacPDstandardNth13gt22,kmadd(gtu12,JacPDstandardNth21gt22,kmadd(gtu22,JacPDstandardNth22gt22,kmadd(gtu23,JacPDstandardNth23gt22,kmadd(gtu13,JacPDstandardNth31gt22,kmadd(gtu33,JacPDstandardNth33gt22,kmul(gtu23,JacPDstandardNth32gt22))))))))),kmadd(Gtl212,Xtn1,kmadd(Gtl222,Xtn2,kmul(Gtl223,Xtn3))))))))))))))))))))));
    
    CCTK_REAL_VEC Rt23 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt112,Gtlu131,kmadd(Gt122,Gtlu132,kmadd(Gt123,Gtlu133,kmadd(Gt113,Gtlu211,kmadd(Gt123,Gtlu212,kmadd(Gt133,Gtlu213,kmadd(Gt213,Gtlu221,kmadd(Gt223,Gtlu222,kmadd(Gt233,Gtlu223,kmadd(Gt212,Gtlu231,kmadd(Gt313,Gtlu231,kmadd(Gt222,Gtlu232,kmadd(Gt323,Gtlu232,kmadd(Gt223,Gtlu233,kmadd(Gt333,Gtlu233,kmadd(Gt112,Gtlu311,kmadd(Gt122,Gtlu312,kmadd(Gt123,Gtlu313,kmadd(Gt212,Gtlu321,kmadd(Gt222,Gtlu322,kmadd(Gt223,Gtlu323,kmadd(ToReal(2),kmul(Gt312,Gtlu331),kmadd(ToReal(2),kmul(Gt322,Gtlu332),kmadd(ToReal(2),kmul(Gt323,Gtlu333),kmadd(ToReal(0.5),kmadd(gt13L,JacPDstandardNth2Xt1,kmadd(gt23L,JacPDstandardNth2Xt2,kmul(gt33L,JacPDstandardNth2Xt3))),knmsub(ToReal(0.5),kmadd(gtu11,JacPDstandardNth11gt23,kmadd(gtu12,JacPDstandardNth12gt23,kmadd(gtu13,JacPDstandardNth13gt23,kmadd(gtu12,JacPDstandardNth21gt23,kmadd(gtu22,JacPDstandardNth22gt23,kmadd(gtu23,JacPDstandardNth23gt23,kmadd(gtu13,JacPDstandardNth31gt23,kmadd(gtu33,JacPDstandardNth33gt23,kmul(gtu23,JacPDstandardNth32gt23))))))))),kmadd(ToReal(0.5),kmadd(gt12L,JacPDstandardNth3Xt1,kmadd(gt22L,JacPDstandardNth3Xt2,kmul(gt23L,JacPDstandardNth3Xt3))),kmadd(ToReal(0.5),kmadd(Gtl213,Xtn1,kmadd(Gtl223,Xtn2,kmul(Gtl233,Xtn3))),kmul(kmadd(Gtl312,Xtn1,kmadd(Gtl322,Xtn2,kmul(Gtl323,Xtn3))),ToReal(0.5))))))))))))))))))))))))))))));
    
    CCTK_REAL_VEC Rt33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt113,Gtlu131,kmadd(Gt123,Gtlu132,kmadd(Gt133,Gtlu133,kmadd(Gt213,Gtlu231,kmadd(Gt223,Gtlu232,kmadd(Gt233,Gtlu233,kmadd(ToReal(2),kmul(Gt113,Gtlu311),kmadd(ToReal(2),kmul(Gt123,Gtlu312),kmadd(ToReal(2),kmul(Gt133,Gtlu313),kmadd(ToReal(2),kmul(Gt213,Gtlu321),kmadd(ToReal(2),kmul(Gt223,Gtlu322),kmadd(ToReal(2),kmul(Gt233,Gtlu323),kmadd(ToReal(3),kmul(Gt313,Gtlu331),kmadd(ToReal(3),kmul(Gt323,Gtlu332),kmadd(ToReal(3),kmul(Gt333,Gtlu333),knmsub(ToReal(0.5),kmadd(gtu11,JacPDstandardNth11gt33,kmadd(gtu12,JacPDstandardNth12gt33,kmadd(gtu13,JacPDstandardNth13gt33,kmadd(gtu12,JacPDstandardNth21gt33,kmadd(gtu22,JacPDstandardNth22gt33,kmadd(gtu23,JacPDstandardNth23gt33,kmadd(gtu13,JacPDstandardNth31gt33,kmadd(gtu33,JacPDstandardNth33gt33,kmul(gtu23,JacPDstandardNth32gt33))))))))),kmadd(gt13L,JacPDstandardNth3Xt1,kmadd(gt23L,JacPDstandardNth3Xt2,kmadd(gt33L,JacPDstandardNth3Xt3,kmadd(Gtl313,Xtn1,kmadd(Gtl323,Xtn2,kmul(Gtl333,Xtn3))))))))))))))))))))));
    
    CCTK_REAL_VEC fac1 CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod == 
      1,kdiv(ToReal(-0.5),phiL),ToReal(1));
    
    CCTK_REAL_VEC cdphi1 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,JacPDstandardNth1phi);
    
    CCTK_REAL_VEC cdphi2 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,JacPDstandardNth2phi);
    
    CCTK_REAL_VEC cdphi3 CCTK_ATTRIBUTE_UNUSED = 
      kmul(fac1,JacPDstandardNth3phi);
    
    CCTK_REAL_VEC fac2 CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod == 
      1,kdiv(ToReal(0.5),kmul(phiL,phiL)),ToReal(0));
    
    CCTK_REAL_VEC cdphi211 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(fac1,ksub(JacPDstandardNth11phi,kmadd(Gt111,JacPDstandardNth1phi,kmadd(Gt311,JacPDstandardNth3phi,kmul(Gt211,JacPDstandardNth2phi)))),kmul(fac2,kmul(JacPDstandardNth1phi,JacPDstandardNth1phi)));
    
    CCTK_REAL_VEC cdphi212 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(fac2,kmul(JacPDstandardNth1phi,JacPDstandardNth2phi),kmul(fac1,ksub(JacPDstandardNth12phi,kmadd(Gt112,JacPDstandardNth1phi,kmadd(Gt312,JacPDstandardNth3phi,kmul(Gt212,JacPDstandardNth2phi))))));
    
    CCTK_REAL_VEC cdphi213 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(fac2,kmul(JacPDstandardNth1phi,JacPDstandardNth3phi),kmul(fac1,ksub(JacPDstandardNth13phi,kmadd(Gt113,JacPDstandardNth1phi,kmadd(Gt313,JacPDstandardNth3phi,kmul(Gt213,JacPDstandardNth2phi))))));
    
    CCTK_REAL_VEC cdphi222 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(fac1,knmsub(Gt122,JacPDstandardNth1phi,ksub(JacPDstandardNth22phi,kmadd(Gt322,JacPDstandardNth3phi,kmul(Gt222,JacPDstandardNth2phi)))),kmul(fac2,kmul(JacPDstandardNth2phi,JacPDstandardNth2phi)));
    
    CCTK_REAL_VEC cdphi223 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(fac2,kmul(JacPDstandardNth2phi,JacPDstandardNth3phi),kmul(fac1,knmsub(Gt123,JacPDstandardNth1phi,ksub(JacPDstandardNth23phi,kmadd(Gt323,JacPDstandardNth3phi,kmul(Gt223,JacPDstandardNth2phi))))));
    
    CCTK_REAL_VEC cdphi233 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(fac1,knmsub(Gt133,JacPDstandardNth1phi,knmsub(Gt233,JacPDstandardNth2phi,knmsub(Gt333,JacPDstandardNth3phi,JacPDstandardNth33phi))),kmul(fac2,kmul(JacPDstandardNth3phi,JacPDstandardNth3phi)));
    
    CCTK_REAL_VEC Rphi11 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),cdphi211,kmadd(ToReal(-2),kmul(gt11L,kmadd(cdphi211,gtu11,kmadd(ToReal(2),kmul(cdphi212,gtu12),kmadd(ToReal(2),kmul(cdphi213,gtu13),kmadd(cdphi222,gtu22,kmadd(ToReal(2),kmul(cdphi223,gtu23),kmul(cdphi233,gtu33))))))),kmadd(ToReal(-4),kmul(gt11L,kmadd(cdphi1,kmadd(cdphi1,gtu11,kmadd(cdphi2,gtu12,kmul(cdphi3,gtu13))),kmadd(cdphi2,kmadd(cdphi1,gtu12,kmadd(cdphi2,gtu22,kmul(cdphi3,gtu23))),kmul(cdphi3,kmadd(cdphi1,gtu13,kmadd(cdphi2,gtu23,kmul(cdphi3,gtu33))))))),kmul(kmul(cdphi1,cdphi1),ToReal(4)))));
    
    CCTK_REAL_VEC Rphi12 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(4),kmul(cdphi1,cdphi2),kmadd(ToReal(-2),cdphi212,kmadd(ToReal(-2),kmul(gt12L,kmadd(cdphi211,gtu11,kmadd(ToReal(2),kmul(cdphi212,gtu12),kmadd(ToReal(2),kmul(cdphi213,gtu13),kmadd(cdphi222,gtu22,kmadd(ToReal(2),kmul(cdphi223,gtu23),kmul(cdphi233,gtu33))))))),kmul(kmul(gt12L,kmadd(cdphi1,kmadd(cdphi1,gtu11,kmadd(cdphi2,gtu12,kmul(cdphi3,gtu13))),kmadd(cdphi2,kmadd(cdphi1,gtu12,kmadd(cdphi2,gtu22,kmul(cdphi3,gtu23))),kmul(cdphi3,kmadd(cdphi1,gtu13,kmadd(cdphi2,gtu23,kmul(cdphi3,gtu33))))))),ToReal(-4)))));
    
    CCTK_REAL_VEC Rphi13 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),cdphi213,kmadd(ToReal(4),kmul(cdphi1,cdphi3),kmadd(ToReal(-2),kmul(gt13L,kmadd(cdphi211,gtu11,kmadd(ToReal(2),kmul(cdphi212,gtu12),kmadd(ToReal(2),kmul(cdphi213,gtu13),kmadd(cdphi222,gtu22,kmadd(ToReal(2),kmul(cdphi223,gtu23),kmul(cdphi233,gtu33))))))),kmul(kmul(gt13L,kmadd(cdphi1,kmadd(cdphi1,gtu11,kmadd(cdphi2,gtu12,kmul(cdphi3,gtu13))),kmadd(cdphi2,kmadd(cdphi1,gtu12,kmadd(cdphi2,gtu22,kmul(cdphi3,gtu23))),kmul(cdphi3,kmadd(cdphi1,gtu13,kmadd(cdphi2,gtu23,kmul(cdphi3,gtu33))))))),ToReal(-4)))));
    
    CCTK_REAL_VEC Rphi22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),cdphi222,kmadd(ToReal(-2),kmul(gt22L,kmadd(cdphi211,gtu11,kmadd(ToReal(2),kmul(cdphi212,gtu12),kmadd(ToReal(2),kmul(cdphi213,gtu13),kmadd(cdphi222,gtu22,kmadd(ToReal(2),kmul(cdphi223,gtu23),kmul(cdphi233,gtu33))))))),kmadd(ToReal(-4),kmul(gt22L,kmadd(cdphi1,kmadd(cdphi1,gtu11,kmadd(cdphi2,gtu12,kmul(cdphi3,gtu13))),kmadd(cdphi2,kmadd(cdphi1,gtu12,kmadd(cdphi2,gtu22,kmul(cdphi3,gtu23))),kmul(cdphi3,kmadd(cdphi1,gtu13,kmadd(cdphi2,gtu23,kmul(cdphi3,gtu33))))))),kmul(kmul(cdphi2,cdphi2),ToReal(4)))));
    
    CCTK_REAL_VEC Rphi23 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),cdphi223,kmadd(ToReal(4),kmul(cdphi2,cdphi3),kmadd(ToReal(-2),kmul(gt23L,kmadd(cdphi211,gtu11,kmadd(ToReal(2),kmul(cdphi212,gtu12),kmadd(ToReal(2),kmul(cdphi213,gtu13),kmadd(cdphi222,gtu22,kmadd(ToReal(2),kmul(cdphi223,gtu23),kmul(cdphi233,gtu33))))))),kmul(kmul(gt23L,kmadd(cdphi1,kmadd(cdphi1,gtu11,kmadd(cdphi2,gtu12,kmul(cdphi3,gtu13))),kmadd(cdphi2,kmadd(cdphi1,gtu12,kmadd(cdphi2,gtu22,kmul(cdphi3,gtu23))),kmul(cdphi3,kmadd(cdphi1,gtu13,kmadd(cdphi2,gtu23,kmul(cdphi3,gtu33))))))),ToReal(-4)))));
    
    CCTK_REAL_VEC Rphi33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(ToReal(-2),cdphi233,kmadd(ToReal(-2),kmul(gt33L,kmadd(cdphi211,gtu11,kmadd(ToReal(2),kmul(cdphi212,gtu12),kmadd(ToReal(2),kmul(cdphi213,gtu13),kmadd(cdphi222,gtu22,kmadd(ToReal(2),kmul(cdphi223,gtu23),kmul(cdphi233,gtu33))))))),kmadd(ToReal(-4),kmul(gt33L,kmadd(cdphi1,kmadd(cdphi1,gtu11,kmadd(cdphi2,gtu12,kmul(cdphi3,gtu13))),kmadd(cdphi2,kmadd(cdphi1,gtu12,kmadd(cdphi2,gtu22,kmul(cdphi3,gtu23))),kmul(cdphi3,kmadd(cdphi1,gtu13,kmadd(cdphi2,gtu23,kmul(cdphi3,gtu33))))))),kmul(kmul(cdphi3,cdphi3),ToReal(4)))));
    
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
    
    CCTK_REAL_VEC R11 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi11,Rt11);
    
    CCTK_REAL_VEC R12 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi12,Rt12);
    
    CCTK_REAL_VEC R13 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi13,Rt13);
    
    CCTK_REAL_VEC R22 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi22,Rt22);
    
    CCTK_REAL_VEC R23 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi23,Rt23);
    
    CCTK_REAL_VEC R33 CCTK_ATTRIBUTE_UNUSED = kadd(Rphi33,Rt33);
    
    CCTK_REAL_VEC trS CCTK_ATTRIBUTE_UNUSED = 
      kmul(em4phi,kmadd(eTxxL,gtu11,kmadd(ToReal(2),kmul(eTxyL,gtu12),kmadd(ToReal(2),kmul(eTxzL,gtu13),kmadd(eTyyL,gtu22,kmadd(ToReal(2),kmul(eTyzL,gtu23),kmul(eTzzL,gtu33)))))));
    
    CCTK_REAL_VEC Ats11 CCTK_ATTRIBUTE_UNUSED = 
      ksub(kmadd(ToReal(4),kmul(cdphi1,JacPDstandardNth1alpha),kmadd(Gt111,JacPDstandardNth1alpha,kmadd(Gt211,JacPDstandardNth2alpha,kmadd(Gt311,JacPDstandardNth3alpha,kmul(alphaL,R11))))),JacPDstandardNth11alpha);
    
    CCTK_REAL_VEC Ats12 CCTK_ATTRIBUTE_UNUSED = 
      ksub(kmadd(Gt112,JacPDstandardNth1alpha,kmadd(Gt212,JacPDstandardNth2alpha,kmadd(ToReal(2),kmadd(cdphi2,JacPDstandardNth1alpha,kmul(cdphi1,JacPDstandardNth2alpha)),kmadd(Gt312,JacPDstandardNth3alpha,kmul(alphaL,R12))))),JacPDstandardNth12alpha);
    
    CCTK_REAL_VEC Ats13 CCTK_ATTRIBUTE_UNUSED = 
      ksub(kmadd(Gt113,JacPDstandardNth1alpha,kmadd(Gt213,JacPDstandardNth2alpha,kmadd(Gt313,JacPDstandardNth3alpha,kmadd(ToReal(2),kmadd(cdphi3,JacPDstandardNth1alpha,kmul(cdphi1,JacPDstandardNth3alpha)),kmul(alphaL,R13))))),JacPDstandardNth13alpha);
    
    CCTK_REAL_VEC Ats22 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt122,JacPDstandardNth1alpha,ksub(kmadd(ToReal(4),kmul(cdphi2,JacPDstandardNth2alpha),kmadd(Gt222,JacPDstandardNth2alpha,kmadd(Gt322,JacPDstandardNth3alpha,kmul(alphaL,R22)))),JacPDstandardNth22alpha));
    
    CCTK_REAL_VEC Ats23 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt123,JacPDstandardNth1alpha,ksub(kmadd(Gt223,JacPDstandardNth2alpha,kmadd(Gt323,JacPDstandardNth3alpha,kmadd(ToReal(2),kmadd(cdphi3,JacPDstandardNth2alpha,kmul(cdphi2,JacPDstandardNth3alpha)),kmul(alphaL,R23)))),JacPDstandardNth23alpha));
    
    CCTK_REAL_VEC Ats33 CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Gt133,JacPDstandardNth1alpha,kmadd(Gt233,JacPDstandardNth2alpha,ksub(kmadd(ToReal(4),kmul(cdphi3,JacPDstandardNth3alpha),kmadd(Gt333,JacPDstandardNth3alpha,kmul(alphaL,R33))),JacPDstandardNth33alpha)));
    
    CCTK_REAL_VEC trAts CCTK_ATTRIBUTE_UNUSED = 
      kmadd(Ats11,gu11,kmadd(ToReal(2),kmul(Ats12,gu12),kmadd(ToReal(2),kmul(Ats13,gu13),kmadd(Ats22,gu22,kmadd(ToReal(2),kmul(Ats23,gu23),kmul(Ats33,gu33))))));
    
    CCTK_REAL_VEC At11rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmadd(At11L,trKL,kmul(kmadd(At11L,Atm11,kmadd(At12L,Atm21,kmul(At13L,Atm31))),ToReal(-2))),kmadd(ToReal(2),kmul(At11L,JacPDstandardNth1beta1),kmadd(ToReal(2),kmul(At12L,JacPDstandardNth1beta2),kmadd(ToReal(2),kmul(At13L,JacPDstandardNth1beta3),kmadd(ToReal(-0.666666666666666666666666666667),kmul(At11L,kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3))),kmadd(em4phi,kmadd(ToReal(-0.333333333333333333333333333333),kmul(g11,trAts),Ats11),kmul(kmul(ToReal(3.14159265358979323846264338328),kmul(alphaL,kmul(em4phi,kmadd(ToReal(-0.333333333333333333333333333333),kmul(g11,trS),eTxxL)))),ToReal(-8))))))));
    
    CCTK_REAL_VEC At12rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmadd(At12L,trKL,kmul(kmadd(At11L,Atm12,kmadd(At12L,Atm22,kmul(At13L,Atm32))),ToReal(-2))),kmadd(At12L,JacPDstandardNth1beta1,kmadd(At22L,JacPDstandardNth1beta2,kmadd(At23L,JacPDstandardNth1beta3,kmadd(At11L,JacPDstandardNth2beta1,kmadd(At12L,JacPDstandardNth2beta2,kmadd(At13L,JacPDstandardNth2beta3,kmadd(ToReal(-0.666666666666666666666666666667),kmul(At12L,kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3))),kmadd(em4phi,kmadd(ToReal(-0.333333333333333333333333333333),kmul(g12,trAts),Ats12),kmul(kmul(ToReal(3.14159265358979323846264338328),kmul(alphaL,kmul(em4phi,kmadd(ToReal(-0.333333333333333333333333333333),kmul(g12,trS),eTxyL)))),ToReal(-8)))))))))));
    
    CCTK_REAL_VEC At13rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmadd(At13L,trKL,kmul(kmadd(At11L,Atm13,kmadd(At12L,Atm23,kmul(At13L,Atm33))),ToReal(-2))),kmadd(At13L,JacPDstandardNth1beta1,kmadd(At23L,JacPDstandardNth1beta2,kmadd(At33L,JacPDstandardNth1beta3,kmadd(At11L,JacPDstandardNth3beta1,kmadd(At12L,JacPDstandardNth3beta2,kmadd(At13L,JacPDstandardNth3beta3,kmadd(ToReal(-0.666666666666666666666666666667),kmul(At13L,kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3))),kmadd(em4phi,kmadd(ToReal(-0.333333333333333333333333333333),kmul(g13,trAts),Ats13),kmul(kmul(ToReal(3.14159265358979323846264338328),kmul(alphaL,kmul(em4phi,kmadd(ToReal(-0.333333333333333333333333333333),kmul(g13,trS),eTxzL)))),ToReal(-8)))))))))));
    
    CCTK_REAL_VEC At22rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmadd(At22L,trKL,kmul(kmadd(At12L,Atm12,kmadd(At22L,Atm22,kmul(At23L,Atm32))),ToReal(-2))),kmadd(ToReal(2),kmul(At12L,JacPDstandardNth2beta1),kmadd(ToReal(2),kmul(At22L,JacPDstandardNth2beta2),kmadd(ToReal(2),kmul(At23L,JacPDstandardNth2beta3),kmadd(ToReal(-0.666666666666666666666666666667),kmul(At22L,kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3))),kmadd(em4phi,kmadd(ToReal(-0.333333333333333333333333333333),kmul(g22,trAts),Ats22),kmul(kmul(ToReal(3.14159265358979323846264338328),kmul(alphaL,kmul(em4phi,kmadd(ToReal(-0.333333333333333333333333333333),kmul(g22,trS),eTyyL)))),ToReal(-8))))))));
    
    CCTK_REAL_VEC At23rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmadd(At23L,trKL,kmul(kmadd(At12L,Atm13,kmadd(At22L,Atm23,kmul(At23L,Atm33))),ToReal(-2))),kmadd(At13L,JacPDstandardNth2beta1,kmadd(At23L,JacPDstandardNth2beta2,kmadd(At33L,JacPDstandardNth2beta3,kmadd(At12L,JacPDstandardNth3beta1,kmadd(At22L,JacPDstandardNth3beta2,kmadd(At23L,JacPDstandardNth3beta3,kmadd(ToReal(-0.666666666666666666666666666667),kmul(At23L,kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3))),kmadd(em4phi,kmadd(ToReal(-0.333333333333333333333333333333),kmul(g23,trAts),Ats23),kmul(kmul(ToReal(3.14159265358979323846264338328),kmul(alphaL,kmul(em4phi,kmadd(ToReal(-0.333333333333333333333333333333),kmul(g23,trS),eTyzL)))),ToReal(-8)))))))))));
    
    CCTK_REAL_VEC At33rhsL CCTK_ATTRIBUTE_UNUSED = 
      kmadd(alphaL,kmadd(At33L,trKL,kmul(kmadd(At13L,Atm13,kmadd(At23L,Atm23,kmul(At33L,Atm33))),ToReal(-2))),kmadd(ToReal(2),kmul(At13L,JacPDstandardNth3beta1),kmadd(ToReal(2),kmul(At23L,JacPDstandardNth3beta2),kmadd(ToReal(2),kmul(At33L,JacPDstandardNth3beta3),kmadd(ToReal(-0.666666666666666666666666666667),kmul(At33L,kadd(JacPDstandardNth1beta1,kadd(JacPDstandardNth2beta2,JacPDstandardNth3beta3))),kmadd(em4phi,kmadd(ToReal(-0.333333333333333333333333333333),kmul(g33,trAts),Ats33),kmul(kmul(ToReal(3.14159265358979323846264338328),kmul(alphaL,kmul(em4phi,kmadd(ToReal(-0.333333333333333333333333333333),kmul(g33,trS),eTzzL)))),ToReal(-8))))))));
    /* Copy local copies back to grid functions */
    vec_store_partial_prepare(i,vecimin,vecimax);
    vec_store_nta_partial(At11rhs[index],At11rhsL);
    vec_store_nta_partial(At12rhs[index],At12rhsL);
    vec_store_nta_partial(At13rhs[index],At13rhsL);
    vec_store_nta_partial(At22rhs[index],At22rhsL);
    vec_store_nta_partial(At23rhs[index],At23rhsL);
    vec_store_nta_partial(At33rhs[index],At33rhsL);
  }
  CCTK_ENDLOOP3STR(ML_BSSN_RHS2);
}
extern "C" void ML_BSSN_RHS2(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_RHS2_Body");
  }
  if (cctk_iteration % ML_BSSN_RHS2_calc_every != ML_BSSN_RHS2_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ML_BSSN::ML_curv",
    "ML_BSSN::ML_curvrhs",
    "ML_BSSN::ML_Gamma",
    "ML_BSSN::ML_lapse",
    "ML_BSSN::ML_log_confac",
    "ML_BSSN::ML_metric",
    "ML_BSSN::ML_shift",
    "ML_BSSN::ML_trace_curv"};
  AssertGroupStorage(cctkGH, "ML_BSSN_RHS2", 8, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_RHS2", 1, 1, 1);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "ML_BSSN_RHS2", 2, 2, 2);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, ML_BSSN_RHS2_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_RHS2_Body");
  }
}

} // namespace ML_BSSN
