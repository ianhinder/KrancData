/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_BSSN_OpenCL_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::At11"),  CCTK_VarIndex("ML_BSSN_OpenCL::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::At12"),  CCTK_VarIndex("ML_BSSN_OpenCL::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::At13"),  CCTK_VarIndex("ML_BSSN_OpenCL::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::At22"),  CCTK_VarIndex("ML_BSSN_OpenCL::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::At23"),  CCTK_VarIndex("ML_BSSN_OpenCL::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::At33"),  CCTK_VarIndex("ML_BSSN_OpenCL::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::A"),  CCTK_VarIndex("ML_BSSN_OpenCL::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::B1"),  CCTK_VarIndex("ML_BSSN_OpenCL::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::B2"),  CCTK_VarIndex("ML_BSSN_OpenCL::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::B3"),  CCTK_VarIndex("ML_BSSN_OpenCL::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::Xt1"),  CCTK_VarIndex("ML_BSSN_OpenCL::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::Xt2"),  CCTK_VarIndex("ML_BSSN_OpenCL::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::Xt3"),  CCTK_VarIndex("ML_BSSN_OpenCL::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::alpha"),  CCTK_VarIndex("ML_BSSN_OpenCL::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::phi"),  CCTK_VarIndex("ML_BSSN_OpenCL::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::gt11"),  CCTK_VarIndex("ML_BSSN_OpenCL::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::gt12"),  CCTK_VarIndex("ML_BSSN_OpenCL::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::gt13"),  CCTK_VarIndex("ML_BSSN_OpenCL::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::gt22"),  CCTK_VarIndex("ML_BSSN_OpenCL::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::gt23"),  CCTK_VarIndex("ML_BSSN_OpenCL::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::gt33"),  CCTK_VarIndex("ML_BSSN_OpenCL::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::beta1"),  CCTK_VarIndex("ML_BSSN_OpenCL::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::beta2"),  CCTK_VarIndex("ML_BSSN_OpenCL::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::beta3"),  CCTK_VarIndex("ML_BSSN_OpenCL::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_OpenCL::trK"),  CCTK_VarIndex("ML_BSSN_OpenCL::trKrhs"));
  
  /* Register all the evolved Array functions with MoL */
  return;
}
