/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

extern "C" void ML_BSSN_DGFE_RegisterSymmetries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* array holding symmetry definitions */
  CCTK_INT sym[3];
  
  /* Register symmetries of grid functions */
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::phi");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::gt11");
  
  sym[0] = -1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::gt12");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::gt13");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::gt22");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::gt23");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::gt33");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::Xt1");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::Xt2");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::Xt3");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::trK");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::At11");
  
  sym[0] = -1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::At12");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::At13");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::At22");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::At23");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::At33");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::alpha");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::A");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::beta1");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::beta2");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::beta3");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::B1");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::B2");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::B3");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::H");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::M1");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::M2");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::M3");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::cS");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::cXt1");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::cXt2");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::cXt3");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "ML_BSSN_DGFE::cA");
  
}
