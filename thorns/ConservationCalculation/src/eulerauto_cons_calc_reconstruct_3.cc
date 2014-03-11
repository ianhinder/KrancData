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

namespace ConservationCalculation {

extern "C" void eulerauto_cons_calc_reconstruct_3_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % eulerauto_cons_calc_reconstruct_3_calc_every != eulerauto_cons_calc_reconstruct_3_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ConservationCalculation::p_lr_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ConservationCalculation::p_lr_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ConservationCalculation::rho_lr_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ConservationCalculation::rho_lr_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ConservationCalculation::v1_lr_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ConservationCalculation::v1_lr_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ConservationCalculation::v2_lr_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ConservationCalculation::v2_lr_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ConservationCalculation::v3_lr_group","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ConservationCalculation::v3_lr_group.");
  return;
}

static void eulerauto_cons_calc_reconstruct_3_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  const CCTK_REAL p1o1 CCTK_ATTRIBUTE_UNUSED = 1;
  const CCTK_REAL p1odx CCTK_ATTRIBUTE_UNUSED = INV(dx);
  const CCTK_REAL p1ody CCTK_ATTRIBUTE_UNUSED = INV(dy);
  const CCTK_REAL p1odz CCTK_ATTRIBUTE_UNUSED = INV(dz);
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
  CCTK_LOOP3(eulerauto_cons_calc_reconstruct_3,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL pL CCTK_ATTRIBUTE_UNUSED = p[index];
    CCTK_REAL rhoL CCTK_ATTRIBUTE_UNUSED = rho[index];
    CCTK_REAL v1L CCTK_ATTRIBUTE_UNUSED = v1[index];
    CCTK_REAL v2L CCTK_ATTRIBUTE_UNUSED = v2[index];
    CCTK_REAL v3L CCTK_ATTRIBUTE_UNUSED = v3[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    const CCTK_REAL DiffPlus3p CCTK_ATTRIBUTE_UNUSED = DiffPlus3(&p[index]);
    const CCTK_REAL DiffMinus3p CCTK_ATTRIBUTE_UNUSED = DiffMinus3(&p[index]);
    const CCTK_REAL DiffPlus3rho CCTK_ATTRIBUTE_UNUSED = DiffPlus3(&rho[index]);
    const CCTK_REAL DiffMinus3rho CCTK_ATTRIBUTE_UNUSED = DiffMinus3(&rho[index]);
    const CCTK_REAL DiffPlus3v1 CCTK_ATTRIBUTE_UNUSED = DiffPlus3(&v1[index]);
    const CCTK_REAL DiffMinus3v1 CCTK_ATTRIBUTE_UNUSED = DiffMinus3(&v1[index]);
    const CCTK_REAL DiffPlus3v2 CCTK_ATTRIBUTE_UNUSED = DiffPlus3(&v2[index]);
    const CCTK_REAL DiffMinus3v2 CCTK_ATTRIBUTE_UNUSED = DiffMinus3(&v2[index]);
    const CCTK_REAL DiffPlus3v3 CCTK_ATTRIBUTE_UNUSED = DiffPlus3(&v3[index]);
    const CCTK_REAL DiffMinus3v3 CCTK_ATTRIBUTE_UNUSED = DiffMinus3(&v3[index]);
    /* Calculate temporaries and grid functions */
    CCTK_REAL slopeL CCTK_ATTRIBUTE_UNUSED = DiffMinus3rho;
    
    CCTK_REAL slopeR CCTK_ATTRIBUTE_UNUSED = DiffPlus3rho;
    
    CCTK_REAL slope CCTK_ATTRIBUTE_UNUSED = VanLeer(slopeL,slopeR);
    
    CCTK_REAL rhoLeftL CCTK_ATTRIBUTE_UNUSED = rhoL - 0.5*slope;
    
    CCTK_REAL rhoRightL CCTK_ATTRIBUTE_UNUSED = rhoL + 0.5*slope;
    
    slopeL = DiffMinus3v1;
    
    slopeR = DiffPlus3v1;
    
    slope = VanLeer(slopeL,slopeR);
    
    CCTK_REAL v1LeftL CCTK_ATTRIBUTE_UNUSED = v1L - 0.5*slope;
    
    CCTK_REAL v1RightL CCTK_ATTRIBUTE_UNUSED = v1L + 0.5*slope;
    
    slopeL = DiffMinus3v2;
    
    slopeR = DiffPlus3v2;
    
    slope = VanLeer(slopeL,slopeR);
    
    CCTK_REAL v2LeftL CCTK_ATTRIBUTE_UNUSED = v2L - 0.5*slope;
    
    CCTK_REAL v2RightL CCTK_ATTRIBUTE_UNUSED = v2L + 0.5*slope;
    
    slopeL = DiffMinus3v3;
    
    slopeR = DiffPlus3v3;
    
    slope = VanLeer(slopeL,slopeR);
    
    CCTK_REAL v3LeftL CCTK_ATTRIBUTE_UNUSED = v3L - 0.5*slope;
    
    CCTK_REAL v3RightL CCTK_ATTRIBUTE_UNUSED = v3L + 0.5*slope;
    
    slopeL = DiffMinus3p;
    
    slopeR = DiffPlus3p;
    
    slope = VanLeer(slopeL,slopeR);
    
    CCTK_REAL pLeftL CCTK_ATTRIBUTE_UNUSED = pL - 0.5*slope;
    
    CCTK_REAL pRightL CCTK_ATTRIBUTE_UNUSED = pL + 0.5*slope;
    /* Copy local copies back to grid functions */
    pLeft[index] = pLeftL;
    pRight[index] = pRightL;
    rhoLeft[index] = rhoLeftL;
    rhoRight[index] = rhoRightL;
    v1Left[index] = v1LeftL;
    v1Right[index] = v1RightL;
    v2Left[index] = v2LeftL;
    v2Right[index] = v2RightL;
    v3Left[index] = v3LeftL;
    v3Right[index] = v3RightL;
  }
  CCTK_ENDLOOP3(eulerauto_cons_calc_reconstruct_3);
}
extern "C" void eulerauto_cons_calc_reconstruct_3(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering eulerauto_cons_calc_reconstruct_3_Body");
  }
  if (cctk_iteration % eulerauto_cons_calc_reconstruct_3_calc_every != eulerauto_cons_calc_reconstruct_3_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ConservationCalculation::p_group",
    "ConservationCalculation::p_lr_group",
    "ConservationCalculation::rho_group",
    "ConservationCalculation::rho_lr_group",
    "ConservationCalculation::v1_lr_group",
    "ConservationCalculation::v2_lr_group",
    "ConservationCalculation::v3_lr_group",
    "ConservationCalculation::v_group"};
  AssertGroupStorage(cctkGH, "eulerauto_cons_calc_reconstruct_3", 8, groups);
  
  EnsureStencilFits(cctkGH, "eulerauto_cons_calc_reconstruct_3", 1, 1, 1);
  
  LoopOverInterior(cctkGH, eulerauto_cons_calc_reconstruct_3_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving eulerauto_cons_calc_reconstruct_3_Body");
  }
}

} // namespace ConservationCalculation
