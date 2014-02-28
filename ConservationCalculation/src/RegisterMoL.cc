/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ConservationCalculation_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ConservationCalculation::Den"),  CCTK_VarIndex("ConservationCalculation::Denrhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ConservationCalculation::En"),  CCTK_VarIndex("ConservationCalculation::Enrhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ConservationCalculation::S1"),  CCTK_VarIndex("ConservationCalculation::S1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ConservationCalculation::S2"),  CCTK_VarIndex("ConservationCalculation::S2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ConservationCalculation::S3"),  CCTK_VarIndex("ConservationCalculation::S3rhs"));
  
  /* Register all the evolved Array functions with MoL */
  return;
}
