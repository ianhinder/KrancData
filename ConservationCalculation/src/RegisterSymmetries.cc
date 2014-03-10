/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

extern "C" void ConservationCalculation_RegisterSymmetries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* array holding symmetry definitions */
  CCTK_INT sym[3];
  
  /* Register symmetries of grid functions */
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::Den");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::S1");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::S2");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::S3");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::En");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::rho");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::v1");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::v2");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::v3");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::p");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::DenFlux");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::DenLeft");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::DenRight");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::EnFlux");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::EnLeft");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::EnRight");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::pLeft");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::pRight");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::rhoLeft");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::rhoRight");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::S1Flux");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::S1Left");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::S1Right");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::S2Flux");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::S2Left");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::S2Right");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::S3Flux");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::S3Left");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::S3Right");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::v1Left");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::v1Right");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::v2Left");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::v2Right");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::v3Left");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ConservationCalculation::v3Right");
  
}
