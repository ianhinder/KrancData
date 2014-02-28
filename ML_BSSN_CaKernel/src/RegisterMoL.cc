/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_BSSN_CaKernel_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::At11"),  CCTK_VarIndex("ML_BSSN_CaKernel::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::At12"),  CCTK_VarIndex("ML_BSSN_CaKernel::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::At13"),  CCTK_VarIndex("ML_BSSN_CaKernel::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::At22"),  CCTK_VarIndex("ML_BSSN_CaKernel::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::At23"),  CCTK_VarIndex("ML_BSSN_CaKernel::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::At33"),  CCTK_VarIndex("ML_BSSN_CaKernel::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::A"),  CCTK_VarIndex("ML_BSSN_CaKernel::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::B1"),  CCTK_VarIndex("ML_BSSN_CaKernel::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::B2"),  CCTK_VarIndex("ML_BSSN_CaKernel::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::B3"),  CCTK_VarIndex("ML_BSSN_CaKernel::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::Xt1"),  CCTK_VarIndex("ML_BSSN_CaKernel::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::Xt2"),  CCTK_VarIndex("ML_BSSN_CaKernel::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::Xt3"),  CCTK_VarIndex("ML_BSSN_CaKernel::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::alpha"),  CCTK_VarIndex("ML_BSSN_CaKernel::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::phi"),  CCTK_VarIndex("ML_BSSN_CaKernel::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::gt11"),  CCTK_VarIndex("ML_BSSN_CaKernel::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::gt12"),  CCTK_VarIndex("ML_BSSN_CaKernel::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::gt13"),  CCTK_VarIndex("ML_BSSN_CaKernel::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::gt22"),  CCTK_VarIndex("ML_BSSN_CaKernel::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::gt23"),  CCTK_VarIndex("ML_BSSN_CaKernel::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::gt33"),  CCTK_VarIndex("ML_BSSN_CaKernel::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::beta1"),  CCTK_VarIndex("ML_BSSN_CaKernel::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::beta2"),  CCTK_VarIndex("ML_BSSN_CaKernel::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::beta3"),  CCTK_VarIndex("ML_BSSN_CaKernel::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_CaKernel::trK"),  CCTK_VarIndex("ML_BSSN_CaKernel::trKrhs"));
  
  /* Register all the evolved Array functions with MoL */
  return;
}
