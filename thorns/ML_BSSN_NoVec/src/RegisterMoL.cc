/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_BSSN_NoVec_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::At11"),  CCTK_VarIndex("ML_BSSN_NoVec::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::At12"),  CCTK_VarIndex("ML_BSSN_NoVec::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::At13"),  CCTK_VarIndex("ML_BSSN_NoVec::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::At22"),  CCTK_VarIndex("ML_BSSN_NoVec::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::At23"),  CCTK_VarIndex("ML_BSSN_NoVec::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::At33"),  CCTK_VarIndex("ML_BSSN_NoVec::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::A"),  CCTK_VarIndex("ML_BSSN_NoVec::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::B1"),  CCTK_VarIndex("ML_BSSN_NoVec::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::B2"),  CCTK_VarIndex("ML_BSSN_NoVec::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::B3"),  CCTK_VarIndex("ML_BSSN_NoVec::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::Xt1"),  CCTK_VarIndex("ML_BSSN_NoVec::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::Xt2"),  CCTK_VarIndex("ML_BSSN_NoVec::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::Xt3"),  CCTK_VarIndex("ML_BSSN_NoVec::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::alpha"),  CCTK_VarIndex("ML_BSSN_NoVec::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::phi"),  CCTK_VarIndex("ML_BSSN_NoVec::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::gt11"),  CCTK_VarIndex("ML_BSSN_NoVec::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::gt12"),  CCTK_VarIndex("ML_BSSN_NoVec::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::gt13"),  CCTK_VarIndex("ML_BSSN_NoVec::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::gt22"),  CCTK_VarIndex("ML_BSSN_NoVec::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::gt23"),  CCTK_VarIndex("ML_BSSN_NoVec::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::gt33"),  CCTK_VarIndex("ML_BSSN_NoVec::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::beta1"),  CCTK_VarIndex("ML_BSSN_NoVec::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::beta2"),  CCTK_VarIndex("ML_BSSN_NoVec::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::beta3"),  CCTK_VarIndex("ML_BSSN_NoVec::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN_NoVec::trK"),  CCTK_VarIndex("ML_BSSN_NoVec::trKrhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
