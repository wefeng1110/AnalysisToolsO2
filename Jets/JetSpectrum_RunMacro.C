// Run inside the ALICE RooUnfold environment:
// root -l Jets/JetSpectrum_RunMacro.C
#if defined(__CLING__)
R__LOAD_LIBRARY(libRooUnfold)
#endif

#include "./JetSpectrum_DrawingFunctions.h"
#include "./JetSpectrum_DrawingFunctions.C"

void SetStyle(Bool_t graypalette=kFALSE);

void JetSpectrum_RunMacro() {
  SetStyle();
  snprintf(optionsAnalysis, sizeof(optionsAnalysis), "%s,%s", unfoldingPrior, unfoldingMethod);

  int iDataset = 0;
  int iRadius = 0;
  int unfoldParameterInput = 10;

  Draw_ResponseMatrices_Fluctuations(iDataset, iRadius);
  Draw_ResponseMatrices_detectorResponse(iDataset, iRadius);
  Draw_ResponseMatrices_DetectorAndFluctuationsCombined(iDataset, iRadius, optionsAnalysis);
  Draw_Pt_spectrum_unfolded_singleDataset(iDataset, iRadius, unfoldParameterInput, optionsAnalysis);
}
