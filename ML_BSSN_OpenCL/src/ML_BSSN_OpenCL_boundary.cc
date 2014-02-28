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
#include "cctk_Loop.h"
#include "loopcontrol.h"
#include "OpenCLRunTime.h"
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

extern "C" void ML_BSSN_OpenCL_boundary_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % ML_BSSN_OpenCL_boundary_calc_every != ML_BSSN_OpenCL_boundary_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_OpenCL::ML_curv","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_OpenCL::ML_curv.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_OpenCL::ML_dtlapse","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_OpenCL::ML_dtlapse.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_OpenCL::ML_dtshift","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_OpenCL::ML_dtshift.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_OpenCL::ML_Gamma","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_OpenCL::ML_Gamma.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_OpenCL::ML_lapse","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_OpenCL::ML_lapse.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_OpenCL::ML_log_confac","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_OpenCL::ML_log_confac.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_OpenCL::ML_metric","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_OpenCL::ML_metric.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_OpenCL::ML_shift","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_OpenCL::ML_shift.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GenericFD_GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN_OpenCL::ML_trace_curv","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for ML_BSSN_OpenCL::ML_trace_curv.");
  return;
}

static void ML_BSSN_OpenCL_boundary_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  const char* const source =
  "\n"
  "/* Include user-supplied include files */\n"
  "\n"
  "/* Initialise finite differencing variables */\n"
  "const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;\n"
  "const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);\n"
  "const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);\n"
  "const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;\n"
  "const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;\n"
  "const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;\n"
  "const CCTK_REAL_VEC dx CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(0));\n"
  "const CCTK_REAL_VEC dy CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(1));\n"
  "const CCTK_REAL_VEC dz CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_SPACE(2));\n"
  "const CCTK_REAL_VEC dt CCTK_ATTRIBUTE_UNUSED = ToReal(CCTK_DELTA_TIME);\n"
  "const CCTK_REAL_VEC t CCTK_ATTRIBUTE_UNUSED = ToReal(cctk_time);\n"
  "const CCTK_REAL_VEC dxi CCTK_ATTRIBUTE_UNUSED = INV(dx);\n"
  "const CCTK_REAL_VEC dyi CCTK_ATTRIBUTE_UNUSED = INV(dy);\n"
  "const CCTK_REAL_VEC dzi CCTK_ATTRIBUTE_UNUSED = INV(dz);\n"
  "const CCTK_REAL_VEC khalf CCTK_ATTRIBUTE_UNUSED = ToReal(0.5);\n"
  "const CCTK_REAL_VEC kthird CCTK_ATTRIBUTE_UNUSED = \n"
  "  ToReal(0.333333333333333333333333333333);\n"
  "const CCTK_REAL_VEC ktwothird CCTK_ATTRIBUTE_UNUSED = \n"
  "  ToReal(0.666666666666666666666666666667);\n"
  "const CCTK_REAL_VEC kfourthird CCTK_ATTRIBUTE_UNUSED = \n"
  "  ToReal(1.33333333333333333333333333333);\n"
  "const CCTK_REAL_VEC hdxi CCTK_ATTRIBUTE_UNUSED = \n"
  "  kmul(dxi,ToReal(0.5));\n"
  "const CCTK_REAL_VEC hdyi CCTK_ATTRIBUTE_UNUSED = \n"
  "  kmul(dyi,ToReal(0.5));\n"
  "const CCTK_REAL_VEC hdzi CCTK_ATTRIBUTE_UNUSED = \n"
  "  kmul(dzi,ToReal(0.5));\n"
  "\n"
  "/* Initialize predefined quantities */\n"
  "const CCTK_REAL_VEC p1o12dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dx);\n"
  "const CCTK_REAL_VEC p1o12dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dy);\n"
  "const CCTK_REAL_VEC p1o12dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0833333333333333333333333333333),dz);\n"
  "const CCTK_REAL_VEC p1o144dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dy,dx));\n"
  "const CCTK_REAL_VEC p1o144dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dz,dx));\n"
  "const CCTK_REAL_VEC p1o144dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.00694444444444444444444444444444),kmul(dz,dy));\n"
  "const CCTK_REAL_VEC p1o24dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dx);\n"
  "const CCTK_REAL_VEC p1o24dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dy);\n"
  "const CCTK_REAL_VEC p1o24dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.0416666666666666666666666666667),dz);\n"
  "const CCTK_REAL_VEC p1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dx);\n"
  "const CCTK_REAL_VEC p1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dy);\n"
  "const CCTK_REAL_VEC p1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.5),dz);\n"
  "const CCTK_REAL_VEC p1o4dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dx);\n"
  "const CCTK_REAL_VEC p1o4dxdy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dy,dx));\n"
  "const CCTK_REAL_VEC p1o4dxdz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dz,dx));\n"
  "const CCTK_REAL_VEC p1o4dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dy);\n"
  "const CCTK_REAL_VEC p1o4dydz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),kmul(dz,dy));\n"
  "const CCTK_REAL_VEC p1o4dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.25),dz);\n"
  "const CCTK_REAL_VEC p1o64dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dx);\n"
  "const CCTK_REAL_VEC p1o64dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dy);\n"
  "const CCTK_REAL_VEC p1o64dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(0.015625),dz);\n"
  "const CCTK_REAL_VEC p1odx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dx);\n"
  "const CCTK_REAL_VEC p1odx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dx,dx));\n"
  "const CCTK_REAL_VEC p1ody CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dy);\n"
  "const CCTK_REAL_VEC p1ody2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dy,dy));\n"
  "const CCTK_REAL_VEC p1odz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dz);\n"
  "const CCTK_REAL_VEC p1odz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),kmul(dz,dz));\n"
  "const CCTK_REAL_VEC pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dx,dx));\n"
  "const CCTK_REAL_VEC pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dy,dy));\n"
  "const CCTK_REAL_VEC pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0833333333333333333333333333333),kmul(dz,dz));\n"
  "const CCTK_REAL_VEC pm1o16dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dx);\n"
  "const CCTK_REAL_VEC pm1o16dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dy);\n"
  "const CCTK_REAL_VEC pm1o16dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.0625),dz);\n"
  "const CCTK_REAL_VEC pm1o2dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dx);\n"
  "const CCTK_REAL_VEC pm1o2dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dy);\n"
  "const CCTK_REAL_VEC pm1o2dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.5),dz);\n"
  "const CCTK_REAL_VEC pm1o4dx CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dx);\n"
  "const CCTK_REAL_VEC pm1o4dy CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dy);\n"
  "const CCTK_REAL_VEC pm1o4dz CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(-0.25),dz);\n"
  "\n"
  "/* Jacobian variable pointers */\n"
  "const bool use_jacobian1 = (!CCTK_IsFunctionAliased(\"MultiPatch_GetMap\") || MultiPatch_GetMap(cctkGH) != jacobian_identity_map)\n"
  "                      && strlen(jacobian_group) > 0;\n"
  "const bool use_jacobian = assume_use_jacobian>=0 ? assume_use_jacobian : use_jacobian1;\n"
  "const bool usejacobian CCTK_ATTRIBUTE_UNUSED = use_jacobian;\n"
  "if (use_jacobian && (strlen(jacobian_derivative_group) == 0))\n"
  "{\n"
  "  CCTK_WARN(1, \"GenericFD::jacobian_group and GenericFD::jacobian_derivative_group must both be set to valid group names\");\n"
  "}\n"
  "\n"
  "const CCTK_REAL* restrict jacobian_ptrs[9];\n"
  "if (use_jacobian) GenericFD_GroupDataPointers(cctkGH, jacobian_group,\n"
  "                                              9, jacobian_ptrs);\n"
  "\n"
  "const CCTK_REAL* restrict const J11 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[0] : 0;\n"
  "const CCTK_REAL* restrict const J12 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[1] : 0;\n"
  "const CCTK_REAL* restrict const J13 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[2] : 0;\n"
  "const CCTK_REAL* restrict const J21 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[3] : 0;\n"
  "const CCTK_REAL* restrict const J22 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[4] : 0;\n"
  "const CCTK_REAL* restrict const J23 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[5] : 0;\n"
  "const CCTK_REAL* restrict const J31 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[6] : 0;\n"
  "const CCTK_REAL* restrict const J32 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[7] : 0;\n"
  "const CCTK_REAL* restrict const J33 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[8] : 0;\n"
  "\n"
  "const CCTK_REAL* restrict jacobian_determinant_ptrs[1] CCTK_ATTRIBUTE_UNUSED;\n"
  "if (use_jacobian && strlen(jacobian_determinant_group) > 0) GenericFD_GroupDataPointers(cctkGH, jacobian_determinant_group,\n"
  "                                              1, jacobian_determinant_ptrs);\n"
  "\n"
  "const CCTK_REAL* restrict const detJ CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_ptrs[0] : 0;\n"
  "\n"
  "const CCTK_REAL* restrict jacobian_inverse_ptrs[9] CCTK_ATTRIBUTE_UNUSED;\n"
  "if (use_jacobian && strlen(jacobian_inverse_group) > 0) GenericFD_GroupDataPointers(cctkGH, jacobian_inverse_group,\n"
  "                                              9, jacobian_inverse_ptrs);\n"
  "\n"
  "const CCTK_REAL* restrict const iJ11 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[0] : 0;\n"
  "const CCTK_REAL* restrict const iJ12 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[1] : 0;\n"
  "const CCTK_REAL* restrict const iJ13 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[2] : 0;\n"
  "const CCTK_REAL* restrict const iJ21 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[3] : 0;\n"
  "const CCTK_REAL* restrict const iJ22 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[4] : 0;\n"
  "const CCTK_REAL* restrict const iJ23 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[5] : 0;\n"
  "const CCTK_REAL* restrict const iJ31 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[6] : 0;\n"
  "const CCTK_REAL* restrict const iJ32 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[7] : 0;\n"
  "const CCTK_REAL* restrict const iJ33 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_inverse_ptrs[8] : 0;\n"
  "\n"
  "const CCTK_REAL* restrict jacobian_derivative_ptrs[18] CCTK_ATTRIBUTE_UNUSED;\n"
  "if (use_jacobian) GenericFD_GroupDataPointers(cctkGH, jacobian_derivative_group,\n"
  "                                              18, jacobian_derivative_ptrs);\n"
  "\n"
  "const CCTK_REAL* restrict const dJ111 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[0] : 0;\n"
  "const CCTK_REAL* restrict const dJ112 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[1] : 0;\n"
  "const CCTK_REAL* restrict const dJ113 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[2] : 0;\n"
  "const CCTK_REAL* restrict const dJ122 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[3] : 0;\n"
  "const CCTK_REAL* restrict const dJ123 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[4] : 0;\n"
  "const CCTK_REAL* restrict const dJ133 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[5] : 0;\n"
  "const CCTK_REAL* restrict const dJ211 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[6] : 0;\n"
  "const CCTK_REAL* restrict const dJ212 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[7] : 0;\n"
  "const CCTK_REAL* restrict const dJ213 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[8] : 0;\n"
  "const CCTK_REAL* restrict const dJ222 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[9] : 0;\n"
  "const CCTK_REAL* restrict const dJ223 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[10] : 0;\n"
  "const CCTK_REAL* restrict const dJ233 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[11] : 0;\n"
  "const CCTK_REAL* restrict const dJ311 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[12] : 0;\n"
  "const CCTK_REAL* restrict const dJ312 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[13] : 0;\n"
  "const CCTK_REAL* restrict const dJ313 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[14] : 0;\n"
  "const CCTK_REAL* restrict const dJ322 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[15] : 0;\n"
  "const CCTK_REAL* restrict const dJ323 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[16] : 0;\n"
  "const CCTK_REAL* restrict const dJ333 CCTK_ATTRIBUTE_UNUSED = use_jacobian ? jacobian_derivative_ptrs[17] : 0;\n"
  "\n"
  "/* Assign local copies of arrays functions */\n"
  "\n"
  "\n"
  "\n"
  "/* Calculate temporaries and arrays functions */\n"
  "\n"
  "/* Copy local copies back to grid functions */\n"
  "\n"
  "/* Loop over the grid points */\n"
  "const int imin0=imin[0];\n"
  "const int imin1=imin[1];\n"
  "const int imin2=imin[2];\n"
  "const int imax0=imax[0];\n"
  "const int imax1=imax[1];\n"
  "const int imax2=imax[2];\n"
  "#pragma omp parallel // reduction(+: vec_iter_counter, vec_op_counter, vec_mem_counter)\n"
  "CCTK_LOOP3STR(ML_BSSN_OpenCL_boundary,\n"
  "  i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,\n"
  "  cctk_ash[0],cctk_ash[1],cctk_ash[2],\n"
  "  vecimin,vecimax, CCTK_REAL_VEC_SIZE)\n"
  "{\n"
  "  const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;\n"
  "  // vec_iter_counter+=CCTK_REAL_VEC_SIZE;\n"
  "  \n"
  "  /* Assign local copies of grid functions */\n"
  "  \n"
  "  \n"
  "  \n"
  "  \n"
  "  /* Include user supplied include files */\n"
  "  \n"
  "  /* Precompute derivatives */\n"
  "  \n"
  "  switch (fdOrder)\n"
  "  {\n"
  "    case 2:\n"
  "    {\n"
  "      break;\n"
  "    }\n"
  "    \n"
  "    case 4:\n"
  "    {\n"
  "      break;\n"
  "    }\n"
  "    default:\n"
  "      CCTK_BUILTIN_UNREACHABLE();\n"
  "  }\n"
  "  \n"
  "  /* Calculate temporaries and grid functions */\n"
  "  CCTK_REAL_VEC phiL CCTK_ATTRIBUTE_UNUSED = IfThen(conformalMethod == \n"
  "    1,ToReal(1),ToReal(0));\n"
  "  \n"
  "  CCTK_REAL_VEC gt11L CCTK_ATTRIBUTE_UNUSED = ToReal(1);\n"
  "  \n"
  "  CCTK_REAL_VEC gt12L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC gt13L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC gt22L CCTK_ATTRIBUTE_UNUSED = ToReal(1);\n"
  "  \n"
  "  CCTK_REAL_VEC gt23L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC gt33L CCTK_ATTRIBUTE_UNUSED = ToReal(1);\n"
  "  \n"
  "  CCTK_REAL_VEC trKL CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC At11L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC At12L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC At13L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC At22L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC At23L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC At33L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC Xt1L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC Xt2L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC Xt3L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC alphaL CCTK_ATTRIBUTE_UNUSED = ToReal(1);\n"
  "  \n"
  "  CCTK_REAL_VEC AL CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC beta1L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC beta2L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC beta3L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC B1L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC B2L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  CCTK_REAL_VEC B3L CCTK_ATTRIBUTE_UNUSED = ToReal(0);\n"
  "  \n"
  "  /* Copy local copies back to grid functions */\n"
  "  vec_store_partial_prepare(i,lc_imin,lc_imax);\n"
  "  vec_store_nta_partial(A[index],AL);\n"
  "  vec_store_nta_partial(alpha[index],alphaL);\n"
  "  vec_store_nta_partial(At11[index],At11L);\n"
  "  vec_store_nta_partial(At12[index],At12L);\n"
  "  vec_store_nta_partial(At13[index],At13L);\n"
  "  vec_store_nta_partial(At22[index],At22L);\n"
  "  vec_store_nta_partial(At23[index],At23L);\n"
  "  vec_store_nta_partial(At33[index],At33L);\n"
  "  vec_store_nta_partial(B1[index],B1L);\n"
  "  vec_store_nta_partial(B2[index],B2L);\n"
  "  vec_store_nta_partial(B3[index],B3L);\n"
  "  vec_store_nta_partial(beta1[index],beta1L);\n"
  "  vec_store_nta_partial(beta2[index],beta2L);\n"
  "  vec_store_nta_partial(beta3[index],beta3L);\n"
  "  vec_store_nta_partial(gt11[index],gt11L);\n"
  "  vec_store_nta_partial(gt12[index],gt12L);\n"
  "  vec_store_nta_partial(gt13[index],gt13L);\n"
  "  vec_store_nta_partial(gt22[index],gt22L);\n"
  "  vec_store_nta_partial(gt23[index],gt23L);\n"
  "  vec_store_nta_partial(gt33[index],gt33L);\n"
  "  vec_store_nta_partial(phi[index],phiL);\n"
  "  vec_store_nta_partial(trK[index],trKL);\n"
  "  vec_store_nta_partial(Xt1[index],Xt1L);\n"
  "  vec_store_nta_partial(Xt2[index],Xt2L);\n"
  "  vec_store_nta_partial(Xt3[index],Xt3L);\n"
  "}\n"
  "CCTK_ENDLOOP3STR(ML_BSSN_OpenCL_boundary);\n"
  ""
  ;
  
  const char* const groups[] = {
    "ML_BSSN_OpenCL::ML_curv",
    "ML_BSSN_OpenCL::ML_dtlapse",
    "ML_BSSN_OpenCL::ML_dtshift",
    "ML_BSSN_OpenCL::ML_Gamma",
    "ML_BSSN_OpenCL::ML_lapse",
    "ML_BSSN_OpenCL::ML_log_confac",
    "ML_BSSN_OpenCL::ML_metric",
    "ML_BSSN_OpenCL::ML_shift",
    "ML_BSSN_OpenCL::ML_trace_curv",
    NULL};
  
  static struct OpenCLKernel *kernel = NULL;
  const char* const sources[] = {differencing, source, NULL};
  OpenCLRunTime_CallKernel(cctkGH, CCTK_THORNSTRING, "ML_BSSN_OpenCL_boundary",
                           sources, groups, NULL, NULL, NULL, -1,
                           imin, imax, &kernel);
  
}

extern "C" void ML_BSSN_OpenCL_boundary(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering ML_BSSN_OpenCL_boundary_Body");
  }
  
  if (cctk_iteration % ML_BSSN_OpenCL_boundary_calc_every != ML_BSSN_OpenCL_boundary_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ML_BSSN_OpenCL::ML_curv",
    "ML_BSSN_OpenCL::ML_dtlapse",
    "ML_BSSN_OpenCL::ML_dtshift",
    "ML_BSSN_OpenCL::ML_Gamma",
    "ML_BSSN_OpenCL::ML_lapse",
    "ML_BSSN_OpenCL::ML_log_confac",
    "ML_BSSN_OpenCL::ML_metric",
    "ML_BSSN_OpenCL::ML_shift",
    "ML_BSSN_OpenCL::ML_trace_curv"};
  GenericFD_AssertGroupStorage(cctkGH, "ML_BSSN_OpenCL_boundary", 9, groups);
  
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
  
  GenericFD_LoopOverBoundaryWithGhosts(cctkGH, ML_BSSN_OpenCL_boundary_Body);
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving ML_BSSN_OpenCL_boundary_Body");
  }
}