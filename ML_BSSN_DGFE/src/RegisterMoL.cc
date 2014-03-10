/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_BSSN_DGFE_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::At11"),  CCTK_VarIndex("ML_BSSN_DGFE::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::At12"),  CCTK_VarIndex("ML_BSSN_DGFE::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::At13"),  CCTK_VarIndex("ML_BSSN_DGFE::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::At22"),  CCTK_VarIndex("ML_BSSN_DGFE::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::At23"),  CCTK_VarIndex("ML_BSSN_DGFE::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::At33"),  CCTK_VarIndex("ML_BSSN_DGFE::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::A"),  CCTK_VarIndex("ML_BSSN_DGFE::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::B1"),  CCTK_VarIndex("ML_BSSN_DGFE::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::B2"),  CCTK_VarIndex("ML_BSSN_DGFE::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::B3"),  CCTK_VarIndex("ML_BSSN_DGFE::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::Xt1"),  CCTK_VarIndex("ML_BSSN_DGFE::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::Xt2"),  CCTK_VarIndex("ML_BSSN_DGFE::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::Xt3"),  CCTK_VarIndex("ML_BSSN_DGFE::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::alpha"),  CCTK_VarIndex("ML_BSSN_DGFE::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::phi"),  CCTK_VarIndex("ML_BSSN_DGFE::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::gt11"),  CCTK_VarIndex("ML_BSSN_DGFE::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::gt12"),  CCTK_VarIndex("ML_BSSN_DGFE::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::gt13"),  CCTK_VarIndex("ML_BSSN_DGFE::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::gt22"),  CCTK_VarIndex("ML_BSSN_DGFE::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::gt23"),  CCTK_VarIndex("ML_BSSN_DGFE::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::gt33"),  CCTK_VarIndex("ML_BSSN_DGFE::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::beta1"),  CCTK_VarIndex("ML_BSSN_DGFE::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::beta2"),  CCTK_VarIndex("ML_BSSN_DGFE::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::beta3"),  CCTK_VarIndex("ML_BSSN_DGFE::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_DGFE::trK"),  CCTK_VarIndex("ML_BSSN_DGFE::trKrhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
