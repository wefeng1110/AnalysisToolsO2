#ifndef JETSPECTRUM_DRAWINGFUNCTIONS_C
#define JETSPECTRUM_DRAWINGFUNCTIONS_C

#include "TStyle.h"
#include "TGraph.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TH1.h"
#include <RooUnfold.h> // one should likely do `aliBuild build RooUnfold` then `alienv enter RooUnfold/latest` as alidist roounfold version is usually quite old
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldSvd.h"
#include "TSVDUnfold.h"

//My Libraries
#include "./JetSpectrum_settings.h"
#include "./JetSpectrum_inputs.h"

#include "./JetSpectrum_ResponseMatrixFunctions.h"
#include "./JetSpectrum_ResponseMatrixFunctions.C"
#include "./JetSpectrum_SpectraGetters.h"
#include "./JetSpectrum_SpectraGetters.C"
#include "./JetSpectrum_Unfolding.h"
#include "./JetSpectrum_Unfolding.C"
#include "./JetSpectrum_EfficiencyPurityGetters.h"
#include "./JetSpectrum_EfficiencyPurityGetters.C"

#include "../Settings/AxisTitles.h"
#include "../Settings/GlobalSettings.h"
#include "../Utilities/AnalysisUtilities.h"
#include "../Utilities/HistogramUtilities.h"
#include "../Utilities/HistogramPlotting.h"
#include "../Utilities/AnalysisUtilities.C" 
#include "../Utilities/HistogramUtilities.C"
#include "../Utilities/HistogramPlotting.C" 
#include "../Utilities/Fits.C"
#include "../Utilities/Fits.h" 

#include<array>
#include <iomanip>
#include <sstream>
#include <string.h>
#include <stdlib.h>     /* abort, NULL */
using namespace std;


/////////////////////////////////////////////////////
/////////////////// Misc utilities //////////////////
/////////////////////////////////////////////////////

void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLineScalePS(1);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
}

void IterationLegend(TString* iterationLegend, int unfoldIterationMin, int unfoldIterationMax, int step){
  const int nUnfoldIteration = std::floor((unfoldIterationMax - unfoldIterationMin)/step) + 1;
  std::stringstream ss;
  ss.precision(2);
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    ss << "k_{unfold} = " << unfoldIterationMax - iUnfoldIteration * step; 
    iterationLegend[iUnfoldIteration] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////// Spectrum plotting functions ///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Draw_Pt_spectrum_raw(int iDataset, int iRadius, std::string options) {

  TH1D* H1D_jetPt_raw;

  if (options.find("noEventNormNorBinWidthScaling") != std::string::npos) {
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_raw, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_bkgCorrected_recBinning(H1D_jetPt_raw, iDataset, iRadius, options);
  }

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_raw");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* yAxisLabel;
  yAxisLabel = texCount;
  if (normaliseDistribsInComparisonPlots) {
    yAxisLabel = texJet_d2Ndptdeta_EventNorm;
  }
  if (options.find("noEventNormNorBinWidthScaling") != std::string::npos) {
    yAxisLabel = texCount;
    *pdfName = *pdfName+(TString)"_noEventNormNorBinWidthScaling";
  }

  Draw_TH1_Histogram(H1D_jetPt_raw, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
}

void Draw_Pt_spectrum_mcp(int iDataset, int iRadius, std::string options) {

  TH1D* H1D_jetPt_mcp_genBinning;
  TH1D* H1D_jetPt_mcp_recBinning;
  TH1D* H1D_jetPt_mcp_collection[2];

  if (options.find("noEventNormNorBinWidthScaling") != std::string::npos) {
    Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp_genBinning, iDataset, iRadius, options);
    Get_Pt_spectrum_mcp_recBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcp_recBinning, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp_genBinning, iDataset, iRadius, options);
    Get_Pt_spectrum_mcp_recBinning(H1D_jetPt_mcp_recBinning, iDataset, iRadius, options);
  }
  
  H1D_jetPt_mcp_collection[0] = H1D_jetPt_mcp_genBinning;
  H1D_jetPt_mcp_collection[1] = H1D_jetPt_mcp_recBinning;


  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_mcp");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* yAxisLabel;
  yAxisLabel = texCount;
  if (normaliseDistribsInComparisonPlots) {
    yAxisLabel = texJet_d2Ndptdeta_EventNorm;
  }
  if (options.find("noEventNormNorBinWidthScaling") != std::string::npos) {
    yAxisLabel = texCount;
    *pdfName = *pdfName+(TString)"_noEventNormNorBinWidthScaling";
  }
  TString genVsRecBinningLegend[2] = {"gen binning", "rec binning"};

  Draw_TH1_Histograms(H1D_jetPt_mcp_collection, genVsRecBinningLegend, 2, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy"); 
}

void Draw_Pt_spectrum_mcdMatched(int iDataset, int iRadius, std::string options) {

  TH1D* H1D_jetPt_mcdMatched;
  if (options.find("noEventNormNorBinWidthScaling") != std::string::npos) {
    Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcdMatched, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_mcdMatched_genBinning(H1D_jetPt_mcdMatched, iDataset, iRadius, options);
  }


  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_mcdMatched");

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* yAxisLabel;
  yAxisLabel = texCount;
  if (normaliseDistribsInComparisonPlots) {
    yAxisLabel = texJet_d2Ndptdeta_EventNorm;
  }
  if (options.find("noEventNormNorBinWidthScaling") != std::string::npos) {
    yAxisLabel = texCount;
    *pdfName = *pdfName+(TString)"_noEventNormNorBinWidthScaling";
  }

  Draw_TH1_Histogram(H1D_jetPt_mcdMatched, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
}


void Draw_Pt_efficiency_jets(int iRadius, std::string options) {
  TH1D* H1D_jetEfficiency[nDatasets];
  bool divideSuccess[nDatasets];
  for (int iDataset = 0; iDataset < nDatasets; ++iDataset) {
    if (useFineBinningTest) {
      divideSuccess[iDataset] = Get_Pt_JetEfficiency_fineBinning(H1D_jetEfficiency[iDataset], iDataset, iRadius, options);
    } else {
      divideSuccess[iDataset] = Get_Pt_JetEfficiency(H1D_jetEfficiency[iDataset], iDataset, iRadius, options);
    }
  }

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_efficiency");
  if (std::all_of(std::begin(divideSuccess), std::end(divideSuccess), [](bool booleanEntry) {return booleanEntry;})){ // checks all entries of divideSuccess are true
    Draw_TH1_Histograms(H1D_jetEfficiency, DatasetsNames, nDatasets, textContext, pdfName, texPtJetGen, texJetEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "efficiency");
  }
}

void Draw_kinematicEfficiency(int iRadius, std::string options) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations[nDatasets];
  TH2D* H2D_jetPtResponseMatrix_detectorResponse[nDatasets];
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning[nDatasets];
  TH1D* H1D_kinematicEfficiency[nDatasets];

  for (int iDataset = 0; iDataset < nDatasets; ++iDataset) {
    TString name_H1D_kinematicEfficiency = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

    Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse[iDataset], iDataset, iRadius, "");

    Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations[iDataset], iDataset, iRadius, "");
    Get_PtResponseMatrix_DetectorAndFluctuationsCombined_preFinalise(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning[iDataset], H2D_jetPtResponseMatrix_detectorResponse[iDataset], H2D_jetPtResponseMatrix_fluctuations[iDataset], iDataset, iRadius, options);

    Get_ResponseMatrix_Pt_KinematicEffiency(H1D_kinematicEfficiency[iDataset], H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning[iDataset], name_H1D_kinematicEfficiency, iRadius);
  }
  TString priorInfo = (TString)unfoldingPrior;

  TString partialUniqueSpecifier = (TString)"R="+Form("%.1f",arrayRadius[iRadius]);
  TString* pdfName = new TString("kinematicEfficiency_"+partialUniqueSpecifier+priorInfo);

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms(H1D_kinematicEfficiency, DatasetsNames, nDatasets, textContext, pdfName, texPtJetGen, texJetKinematicEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
}


void Draw_FakeRatio(int iRadius, std::string options) {
  TH1D* H1D_fakeRatio[nDatasets];

  TString partialUniqueSpecifier = (TString)" R="+Form("%.1f",arrayRadius[iRadius]);

  for (int iDataset = 0; iDataset < nDatasets; ++iDataset) {
    Get_Pt_JetFakes(H1D_fakeRatio[iDataset], iDataset, iRadius, options);
  }
  TString* pdfName = new TString("fakeRatio_"+partialUniqueSpecifier);

  TString textContext(contextCustomOneField(*texDatasetsComparisonCommonDenominator, ""));

  Draw_TH1_Histograms(H1D_fakeRatio, DatasetsNames, nDatasets, textContext, pdfName, texPtJetRec, texFakeRatio, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "");
}

void Draw_ResponseMatrices_Fluctuations(int iDataset, int iRadius) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations;

  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, "");

  TString priorInfo = (TString)unfoldingPrior;

  std::error_code errPDF, errPNG, errEPS;
  CreateDirectoryRecursive((std::string)"pdfFolder/ResponseMatrices", errPDF);
  CreateDirectoryRecursive((std::string)"pngFolder/ResponseMatrices", errPNG);
  CreateDirectoryRecursive((std::string)"epsFolder/ResponseMatrices", errEPS);
  // struct stat st1{};
  // if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
  //     mkdir("pdfFolder/ResponseMatrices", 0700);
  // }
  // struct stat st2{};
  // if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
  //     mkdir("pngFolder/ResponseMatrices", 0700);
  // }
  // struct stat st3{};
  // if (stat("epsFolder/ResponseMatrices", &st3) == -1) {
  //     mkdir("epsFolder/ResponseMatrices", 0700);
  // }

  TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_fluctuationsBackground_"+(TString)"_R="+arrayRadiusPdfName[iRadius]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");
  // TString* pdfNameFullRes_logz = new TString("ResponseMatrices/responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"FullRes_logz");
  TString* pdfName = new TString("ResponseMatrices/responseMatrix_fluctuationsBackground_"+(TString)"_R="+arrayRadiusPdfName[iRadius]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
  // TString* pdfNameFullRes = new TString("ResponseMatrices/responseMatrix_fluctuationsBackground_"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_FullRes");


  TString texCombinedMatrix = contextCustomOneField((TString)"ALICE Performance", ""); // Response matrix - "+(TString)*texEnergy
  TString textContextMatrixDetails = contextCustomFiveFields((TString)"Bkg. fluctuation response ", "", (TString)*texCollisionDataType, (TString)*texEnergyPbPb, contextJetRadius(arrayRadius[iRadius]), "");

  // the matrix natural visualisation is actually the NON transposed histograms, rotated by 90° anti trigonometrically
  TH2D* MatrixResponse;
  TString* xLabel;
  TString* yLabel;
  if (transposeResponseHistogramsInDrawing) {
    MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_fluctuations).Clone("Draw_ResponseMatrices_Fluctuations"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    xLabel = texPtJetFluctCorrectedX;
    yLabel = texPtJetRec;
  } else {
    MatrixResponse = (TH2D*)H2D_jetPtResponseMatrix_fluctuations->Clone("Draw_ResponseMatrices_Fluctuations"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    xLabel = texPtJetRec;
    yLabel = texPtJetFluctCorrectedX;
  }

  double th2ContourCustom[1] = {0.000001}; // hardcoded at 10-6 for now
  int contourNumberCustom = 1;

  // std::array<std::array<float, 2>, 3> drawnWindowYaxianRequest = {{{-999, -999}, {-999, -999}, {1e-6, 1e-1}}}; // {{xmin, xmax}, {ymin, ymax}, {zmin, zmax}} /// put AUTO again after perf figure is done

  // Draw_TH2_Histogram(H2D_jetPtResponseMatrix_fluctuations, textContext, pdfName_logz, texPtJetRec, texPtJetFluctCorrectedX, texCollisionDataInfo, drawnWindow2DAuto, th2ContourCustom, contourNumberCustom, "logz");
  Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_logz, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
}

void Draw_ResponseMatrices_detectorResponse(int iDataset, int iRadius) {

  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  cout << "Draw_ResponseMatrices_detectorResponse 1" << endl;
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius, "");
  cout << "Draw_ResponseMatrices_detectorResponse 2" << endl;

  TString priorInfo = (TString)unfoldingPrior;


  std::error_code errPDF, errPNG, errEPS;
  CreateDirectoryRecursive((std::string)"pdfFolder/ResponseMatrices", errPDF);
  CreateDirectoryRecursive((std::string)"pngFolder/ResponseMatrices", errPNG);
  CreateDirectoryRecursive((std::string)"epsFolder/ResponseMatrices", errEPS);
  // struct stat st1{};
  // if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
  //     mkdir("pdfFolder/ResponseMatrices", 0700);
  // }
  // struct stat st2{};
  // if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
  //     mkdir("pngFolder/ResponseMatrices", 0700);
  // }
  // struct stat st3{};
  // if (stat("epsFolder/ResponseMatrices", &st3) == -1) {
  //     mkdir("epsFolder/ResponseMatrices", 0700);
  // }

  TString* pdfName = new TString("ResponseMatrices/responseMatrix_detectorEffects_"+jetType[iJetType]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
  TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_detectorEffects_"+(TString)"_R="+arrayRadiusPdfName[iRadius]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");

  TString texCombinedMatrix = contextCustomOneField((TString)"ALICE Simulation", ""); // Response matrix - "+(TString)*texEnergy
  TString textContextMatrixDetails = contextCustomFiveFields((TString)"Detector response ", "", (TString)*texCollisionMCType, (TString)*texEnergy, (TString)contextJetRadius(arrayRadius[iRadius]), "");


  // the matrix natural visualisation is actually the NON transposed histograms, rotated by 90° anti trigonometrically
  TH2D* MatrixResponse;
  TString* xLabel;
  TString* yLabel;
  if (transposeResponseHistogramsInDrawing) {
    MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_detectorResponse).Clone("Draw_ResponseMatrices_detectorResponse"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    xLabel = texPtJetGen;
    yLabel = texPtJetRec;
  } else {    MatrixResponse = (TH2D*)H2D_jetPtResponseMatrix_detectorResponse->Clone("Draw_ResponseMatrices_detectorResponse"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    xLabel = texPtJetRec;
    yLabel = texPtJetGen;
  }
  // std::array<std::array<float, 2>, 3> drawnWindowRaymondRequest = {{{-999, -999}, {-999, -999}, {1e-5, 6e-1}}}; // {{xmin, xmax}, {ymin, ymax}, {zmin, zmax}} /// put AUTO again after perf figure is done

  Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
  Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_logz, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
}

void Draw_ResponseMatrices_DetectorAndFluctuationsCombined(int iDataset, int iRadius, std::string options) {

  TH2D* H2D_jetPtResponseMatrix_fluctuations;
  TH2D* H2D_jetPtResponseMatrix_detectorResponse;
  TH2D* H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined;


  Get_PtResponseMatrix_Fluctuations(H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, "");
  Get_PtResponseMatrix_detectorResponse(H2D_jetPtResponseMatrix_detectorResponse, iDataset, iRadius, "");
  Get_PtResponseMatrix_DetectorAndFluctuationsCombined_postFinalise(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, H2D_jetPtResponseMatrix_detectorResponse, H2D_jetPtResponseMatrix_fluctuations, iDataset, iRadius, options);
  // FinaliseResponseMatrix_priorAndNormYslicesAndMergeBins(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, iDataset, iRadius, options);

  TString priorInfo = (TString)unfoldingPrior;


  std::error_code errPDF, errPNG, errEPS;
  CreateDirectoryRecursive((std::string)"pdfFolder/ResponseMatrices", errPDF);
  CreateDirectoryRecursive((std::string)"pngFolder/ResponseMatrices", errPNG);
  CreateDirectoryRecursive((std::string)"epsFolder/ResponseMatrices", errEPS);
  // struct stat st1{};
  // if (stat("pdfFolder/ResponseMatrices", &st1) == -1) {
  //     mkdir("pdfFolder/ResponseMatrices", 0700);
  // }
  // struct stat st2{};
  // if (stat("pngFolder/ResponseMatrices", &st2) == -1) {
  //     mkdir("pngFolder/ResponseMatrices", 0700);
  // }
  // struct stat st3{};
  // if (stat("epsFolder/ResponseMatrices", &st3) == -1) {
  //     mkdir("epsFolder/ResponseMatrices", 0700);
  // }

  TString* pdfName = new TString("ResponseMatrices/responseMatrix_combined"+(TString)"_R="+arrayRadiusPdfName[iRadius]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
  TString* pdfName_logz = new TString("ResponseMatrices/responseMatrix_combined"+(TString)"_R="+arrayRadiusPdfName[iRadius]+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo+"_logz");

  TString texCombinedMatrix = contextCustomOneField((TString)"Combined matrix - "+(TString)*texEnergy, "");
  TString textContextMatrixDetails = contextCustomFourFields((TString)"Detector response: "+(TString)*texCollisionMCType, "", (TString)"Fluctuations response: "+*texCollisionDataType, contextJetRadius(arrayRadius[iRadius]), "");

  // the matrix natural visualisation is actually the NON transposed histograms, rotated by 90° anti trigonometrically
  TH2D* MatrixResponse;
  TString* xLabel;
  TString* yLabel;
  if (transposeResponseHistogramsInDrawing) {
    MatrixResponse = (TH2D*)GetTransposeHistogram(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined).Clone("Draw_ResponseMatrices_DetectorAndFluctuationsCombined"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    xLabel = texPtJetGen;
    yLabel = texPtJetRec;
  } else {
    MatrixResponse = (TH2D*)H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined->Clone("Draw_ResponseMatrices_DetectorAndFluctuationsCombined"+(TString)"_R="+Form("%.1f",arrayRadius[iRadius])+"_"+Datasets[iDataset]+DatasetsNames[iDataset]+"_"+priorInfo);
    xLabel = texPtJetRec;
    yLabel = texPtJetGen;
  }

  Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "");
  Draw_TH2_Histogram(MatrixResponse, textContextMatrixDetails, pdfName_logz, xLabel, yLabel, &texCombinedMatrix, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");

  // save matrix in .root file
  TString fileName = "jetPtResponseMatrixCombined";
  Save_PtResponseMatrix(H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, fileName);
}

void Draw_Pt_spectrum_unfolded_singleDataset(int iDataset, int iRadius, int unfoldParameterInput, std::string options) {
  bool splitTestControlMC = true;

  TH1D* H1D_jetPt_measured;
  TH1D* H1D_jetPt_measured_genBinning;
  TH1D* H1D_jetPt_unfolded;
  TH1D* H1D_jetPt_unfoldedThenRefolded;
  TH1D* H1D_jetPt_unfoldedThenRefolded_RooUnfoldMethod;
  TH1D* H1D_jetPt_mcpFolded;
  TH1D* H1D_jetPt_mcpFolded2;
  TH1D* H1D_jetPt_mcpFoldedThenUnfolded;
  TH1D* H1D_jetPt_unfolded_mcpComp[2];
  TH1D* H1D_jetPt_unfolded_run2Comp_fitRebin[3];
  TH1D* H1D_jetPt_unfolded_run2Comp_shapeComp[2];
  TH1D* H1D_jetPt_unfolded_run2Comp[3];
  TH1D* H1D_jetPt_unfolded_measuredComp[2];
  TH1D* H1D_jetPt_unfolded_refoldedComp[3];
  TH1D* H1D_jetPt_unfolded_mcpFoldedComp[2];
  TH1D* H1D_jetPt_unfolded_mcpFoldedUnfoldedComp[2];
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_mcp_recBinControl;
  TH1D* H1D_jetPt_run2_HannaBossiLauraFile;
  TGraph* Graph_jetPt_run2_MLPaperFile;
  TH1D* H1D_jetPt_run2_MLPaperFile = new TH1D("H1D_jetPt_run2_MLPaperFile", "H1D_jetPt_run2_MLPaperFile", nBinPtJetsGen_run2[iRadius], ptBinsJetsGen_run2[iRadius]);
  TH1D* H1D_jetPt_run2_MLPaperFile_rebinned;
  std::vector<TGraphErrors*> TGraph_jetPt_run2_MLPaperFile_fit = {};
  std::vector<TGraphErrors*> TGraph_jetPt_unfolded_run2Comp_fits = {};
  std::vector<TGraphErrors*> TGraph_xT_run2Comp_fits = {};
  
  TH1D* H1D_jetPt_ratio_mcp;
  TH1D* H1D_jetPt_ratio_run2_fitRebin[2];
  TH1D* H1D_jetPt_ratio_run2_shapeComp[2];
  TH1D* H1D_jetPt_ratio_run2[2];
  TH1D* H1D_jetPt_unfolded_run2Comp_xT[2];
  TH1D* H1D_jetPt_ratio_run2Comp_xT;
  TH1D* H1D_jetPt_unfolded_run2Comp_fits[2];
  TH1D* H1D_jetPt_ratio_run2Comp_fits;
  TH1D* H1D_xT_ratio_run2Comp_fits;
  TH1D* H1D_jetPt_ratio_measured;
  TH1D* H1D_jetPt_ratio_measuredRefolded[2];
  TH1D* H1D_jetPt_ratio_mcpFoldedMcp;
  TH1D* H1D_jetPt_ratio_mcpFoldedUnfoldedMcp;

  TH1D* measuredInput_mcSplitInput;
  TH1D* H1D_jetPt_mcp_mcSplitInput;
  TH1D* H1D_jetPt_mcp_mcSplitInput_recBinning;
  TH1D* H1D_jetPt_unfolded_inputSplitClosure;
  TH1D* H1D_jetPt_unfolded_mcdSplitClosure[2];
  TH1D* H1D_jetPt_ratio_mcdSplitClosure;
  // RUN 2 settings
  if (comparePbPbWithRun2) {
    H1D_jetPt_run2_HannaBossiLauraFile = (TH1D*)((TH1D*)file_O2Analysis_run2ComparisonFileHannaBossiLaura->Get("Bayesian_Unfoldediter15"))->Clone("H1D_jetPt_run2_HannaBossiLauraFile");
    int NcollRun2 = 4619963; // central (see Laura discussion mattermost) 
    H1D_jetPt_run2_HannaBossiLauraFile->Scale(1./NcollRun2);

    double Ncoll;
    if (centralityRange[0] == 00 && centralityRange[1] == 10) {
      // Ncoll = (1780.9+1387.0)/2; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
      Ncoll = (1956+1722+1521+1346)/4; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
    } else if (centralityRange[0] == 50 && centralityRange[1] == 70) {
      // Ncoll = (103.7+46.1)/2; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
      Ncoll = (89.8+39.8)/2; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
    } else {
      cout << "comparison with run2: Ncoll hasn't been calculated for this centrality interval" << endl;
    }
    double sigmaNN = 67.6; // value for sqrt(s) = 5.02 TeV https://arxiv.org/abs/1710.07098
    double T_AA = Ncoll / sigmaNN;
    Graph_jetPt_run2_MLPaperFile = ((TGraph*)((TDirectoryFile*)file_O2Analysis_run2ComparisonFileMLPaper->Get("Figure 3a top R020"))->FindObjectAny("Graph1D_y1")); // https://doi.org/10.1016/j.physletb.2023.138412
    // H1D_jetPt_run2_MLPaperFile = (TH1D*)((TH1D*)(file_O2Analysis_run2ComparisonFileMLPaper->Get("Figure 3a top R020"))->FindObject("Graph1D_y1"))->Clone("H1D_jetPt_run2_MLPaperFile");
    int Ngraph = Graph_jetPt_run2_MLPaperFile->GetN();
    for (int i=0; i < Ngraph; ++i) // setting bin contents to the TGraph values
    {
      double x,y;
      Graph_jetPt_run2_MLPaperFile->GetPoint(i, x, y);
      H1D_jetPt_run2_MLPaperFile->Fill(x, y); // uncertainties are of course screwed up
      int iHist = H1D_jetPt_run2_MLPaperFile->GetXaxis()->FindBin(x);
      H1D_jetPt_run2_MLPaperFile->SetBinError(iHist, Graph_jetPt_run2_MLPaperFile->GetErrorY(i));
    }
    H1D_jetPt_run2_MLPaperFile->Scale(T_AA);
    // now the spectre from the file is 1/N d2N/dpTdeta, instead of 1/T_AA 1/N d2N/dpTdeta
  }

  bool divideSuccessMcp;
  bool divideSuccessRun2_fitRebin[2];
  bool divideSuccessRun2_shapeComp;
  bool divideSuccessRun2[2];
  bool divideSuccessRun2_xt;
  bool divideSuccessRun2_fits; 
  bool divideSuccessRun2_xT_fits; 
  bool divideSuccessMeasured;
  bool divideSuccessMeasuredRefolded;
  bool divideSuccessMcpFoldedMcp;
  bool divideSuccessMcpFoldedUnfoldedMcp;
  bool divideSuccessMcdSplitClosure;
  TString partialUniqueSpecifier;

  if (!useFineBinningTest) {
    Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp, iDataset, iRadius, options);
    if (doClosure_splitMC_mcdUnfoldedVsGen) {
      Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp_mcSplitInput, iDataset, iRadius, options, splitTestControlMC);
    }
  } else {
    Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp, iDataset, iRadius, options);
    if (doClosure_splitMC_mcdUnfoldedVsGen) {
      Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp_mcSplitInput, iDataset, iRadius, options, splitTestControlMC);
    }
  }
  // H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  int unfoldParameter, unfoldParameter_mcSplitClosure;

  partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  TH1D* measuredInput;
  if (!normGenAndMeasByNEvtsForUnfoldingInput) {
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options);
      if (doClosure_splitMC_mcdUnfoldedVsGen) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput_mcSplitInput, iDataset, iRadius, options, splitTestControlMC);
      }
    } else {
      Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options); 
      if (doClosure_splitMC_mcdUnfoldedVsGen) {
        Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput_mcSplitInput, iDataset, iRadius, options, splitTestControlMC);
      }
    }
  } else{
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
      if (doClosure_splitMC_mcdUnfoldedVsGen) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput_mcSplitInput, iDataset, iRadius, options, splitTestControlMC);
      }
    } else {
      Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
      if (doClosure_splitMC_mcdUnfoldedVsGen) {
        Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput_mcSplitInput, iDataset, iRadius, options, splitTestControlMC);
      }
    }
  }

  unfoldParameter = Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options).first;
  // TH1D* H1D_jetPt_unfolded2 = (TH1D*)H1D_jetPt_unfolded->Clone(H1D_jetPt_unfolded->GetName()+(TString)"H1D_jetPt_unfolded2");

  cout << "comparison with measured" << endl; 
  if (!useFineBinningTest) {
    Get_Pt_spectrum_bkgCorrected_genBinning(H1D_jetPt_measured_genBinning, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_bkgCorrected_fineBinning(H1D_jetPt_measured_genBinning, iDataset, iRadius, options);
  }
  H1D_jetPt_unfolded_measuredComp[0] = (TH1D*)H1D_jetPt_measured_genBinning->Clone("H1D_jetPt_measured_genBinning_measuredComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_measuredComp[1] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_measuredComp"+partialUniqueSpecifier);
  H1D_jetPt_ratio_measured = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_ratio_measured"+partialUniqueSpecifier);
  // divideSuccessMeasured = H1D_jetPt_ratio_measured->Divide(H1D_jetPt_measured_genBinning);
  divideSuccessMeasured = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_measured, H1D_jetPt_measured_genBinning);
  

  cout << "comparison with mcp truth" << endl; 
  H1D_jetPt_unfolded_mcpComp[0] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_mcp_mcpComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_mcpComp[1] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
  H1D_jetPt_ratio_mcp = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_ratio_mcp"+partialUniqueSpecifier);
  divideSuccessMcp = H1D_jetPt_ratio_mcp->Divide(H1D_jetPt_mcp);

  cout << "comparison with run2" << endl; 
  std::vector<double> xtBinningVectorRun2 = {};
  std::vector<double> xtBinningVectorRun3 = {};
  if (comparePbPbWithRun2) {
    // comparison with run2 results rebinned using a fit (errors look way underestimated; tsallis function not great aboe 100+ GeV; where should one eval the function inside a bin? probably not just the center)
    H1D_jetPt_unfolded_run2Comp_fitRebin[0] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_run2Comp_fitRebin"+partialUniqueSpecifier);
    double fitPtRange[2] = {ptBinsJetsGen_run2[iRadius][0], ptMaxFit}; //-1 because Tsallis shape only accurate until 120GeV or so
    std::tuple<TF1*, TMatrixDSym, TFitResultPtr> exponentialFitWithLogTransfoResult_run2 = ExponentialFitWithLogTransfo(H1D_jetPt_run2_MLPaperFile, fitPtRange);
    TString histName_run2 = "rebinWithFit_run2"+partialUniqueSpecifier;
    std::pair<TH1D*, TGraphErrors*> pairResult_run2 = RebinWithFit(H1D_jetPt_run2_MLPaperFile, nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius], fitPtRange, histName_run2, exponentialFitWithLogTransfoResult_run2);
    H1D_jetPt_run2_MLPaperFile_rebinned = pairResult_run2.first;
    TGraph_jetPt_run2_MLPaperFile_fit.push_back(pairResult_run2.second);
    H1D_jetPt_unfolded_run2Comp_fitRebin[1] = (TH1D*)H1D_jetPt_run2_MLPaperFile_rebinned->Clone("H1D_jetPt_unfolded_run2_rebinned_fitRebin"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_run2Comp_fitRebin[2] = (TH1D*)H1D_jetPt_run2_MLPaperFile->Clone("H1D_jetPt_unfolded_run2_fitRebin"+partialUniqueSpecifier);
    // H1D_jetPt_unfolded_run2Comp[2] = (TH1D*)H1D_jetPt_run2_HannaBossiLauraFile->Clone("H1D_jetPt_unfolded_run2Comp_HannaBossiLauraFile"+partialUniqueSpecifier);
    H1D_jetPt_ratio_run2_fitRebin[0] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_ratio_run2_fitRebin"+partialUniqueSpecifier);
    // H1D_jetPt_ratio_run2[1] = (TH1D*)H1D_jetPt_run2_HannaBossiLauraFile->Clone("H1D_jetPt_ratio_run2_HannaBossiLauraFile"+partialUniqueSpecifier);
    divideSuccessRun2_fitRebin[0] = H1D_jetPt_ratio_run2_fitRebin[0]->Divide(H1D_jetPt_run2_MLPaperFile_rebinned);
    // divideSuccessRun2[1] = H1D_jetPt_ratio_run2[1]->Divide(H1D_jetPt_unfolded);


    // comparison with run2 results by rebinning Run3 into Run2 bins
    H1D_jetPt_unfolded_run2Comp[0] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_run2Comp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_run2Comp[1] = (TH1D*)H1D_jetPt_unfolded->Rebin(nBinPtJetsGen_run2[iRadius],"H1D_jetPt_unfolded_run2Comp_run2Rebin"+partialUniqueSpecifier, ptBinsJetsGen_run2[iRadius]);
    
    int scalingFactorRebin[9] = {2, 2, 2, 2, 2, 3, 3, 1, 1}; //width of run 2 histogram
    for (auto i = 1; i <= H1D_jetPt_unfolded_run2Comp[1]->GetNbinsX(); i++) {
      H1D_jetPt_unfolded_run2Comp[1]->SetBinContent(i, 1./scalingFactorRebin[i-1]*H1D_jetPt_unfolded_run2Comp[1]->GetBinContent(i));
      H1D_jetPt_unfolded_run2Comp[1]->SetBinError(i, 1./scalingFactorRebin[i-1]*H1D_jetPt_unfolded_run2Comp[1]->GetBinError(i));
    }
    H1D_jetPt_unfolded_run2Comp[2] = (TH1D*)H1D_jetPt_run2_MLPaperFile->Clone("H1D_jetPt_unfolded_run2"+partialUniqueSpecifier);
    // H1D_jetPt_unfolded_run2Comp[2] = (TH1D*)H1D_jetPt_run2_HannaBossiLauraFile->Clone("H1D_jetPt_unfolded_run2Comp_HannaBossiLauraFile"+partialUniqueSpecifier);
    H1D_jetPt_ratio_run2[0] = (TH1D*)H1D_jetPt_unfolded_run2Comp[1]->Clone("H1D_jetPt_ratio_run2"+partialUniqueSpecifier);
    // H1D_jetPt_ratio_run2[1] = (TH1D*)H1D_jetPt_run2_HannaBossiLauraFile->Clone("H1D_jetPt_ratio_run2_HannaBossiLauraFile"+partialUniqueSpecifier);
    divideSuccessRun2[0] = H1D_jetPt_ratio_run2[0]->Divide(H1D_jetPt_run2_MLPaperFile);
    // divideSuccessRun2[1] = H1D_jetPt_ratio_run2[1]->Divide(H1D_jetPt_unfolded);

    H1D_jetPt_unfolded_run2Comp_shapeComp[0] = (TH1D*)H1D_jetPt_unfolded_run2Comp[1]->Clone("H1D_jetPt_unfolded_run2Comp_run3rebinned_shapeComp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_run2Comp_shapeComp[1] = (TH1D*)H1D_jetPt_run2_MLPaperFile->Clone("H1D_jetPt_unfolded_run2Comp_run2rescaled_shapeComp"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_run2Comp_shapeComp[1]->Scale(H1D_jetPt_unfolded_run2Comp_shapeComp[0]->GetBinContent(1)/H1D_jetPt_unfolded_run2Comp_shapeComp[1]->GetBinContent(1));

    H1D_jetPt_ratio_run2_shapeComp[0] = (TH1D*)H1D_jetPt_unfolded_run2Comp_shapeComp[0]->Clone("H1D_jetPt_unfolded_run2Comp_run2rescaled_shapeComp_ratio"+partialUniqueSpecifier);
    divideSuccessRun2_shapeComp = H1D_jetPt_ratio_run2_shapeComp[0]->Divide(H1D_jetPt_unfolded_run2Comp_shapeComp[1]);

    ///////////////////////////////////////
    // xT comparison: xT=2pT/sqrt(s) //////
    ///////////////////////////////////////

    double sqrtS_run2 = 5020; // same unit as pT
    double sqrtS_run3 = 5360; // same unit as pT

    // H1D_jetPt_unfolded_run2Comp_xT[0] = H1D_jetPt_run2_MLPaperFile->Clone("H1D_jetPt_unfolded_run2Comp_xT_run2"+partialUniqueSpecifier);
    // H1D_jetPt_unfolded_run2Comp_xT[1] = H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_run2Comp_xT_run3"+partialUniqueSpecifier);

    for (int iBin = 1; iBin <= H1D_jetPt_run2_MLPaperFile->GetNbinsX()+1; iBin++) {
      xtBinningVectorRun2.push_back(2./sqrtS_run2*H1D_jetPt_run2_MLPaperFile->GetXaxis()->GetBinLowEdge(iBin));
    }
    double* xtBinningRun2 = &xtBinningVectorRun2[0];
    for (int iBin = 1; iBin <= H1D_jetPt_unfolded->GetNbinsX()+1; iBin++) {
      xtBinningVectorRun3.push_back(2./sqrtS_run3*H1D_jetPt_unfolded->GetXaxis()->GetBinLowEdge(iBin));
    }
    double* xtBinningRun3 = &xtBinningVectorRun3[0];

    H1D_jetPt_unfolded_run2Comp_xT[0] = new TH1D("H1D_jetPt_unfolded_run2Comp_xT_run2", "H1D_jetPt_unfolded_run2Comp_xT_run2", H1D_jetPt_run2_MLPaperFile->GetNbinsX(), xtBinningRun2);
    H1D_jetPt_unfolded_run2Comp_xT[1] = new TH1D("H1D_jetPt_unfolded_run2Comp_xT_run3", "H1D_jetPt_unfolded_run2Comp_xT_run3", H1D_jetPt_unfolded->GetNbinsX(), xtBinningRun3);

    double dpt_dxt_run2=sqrtS_run2/2.;
    double dpt_dxt_run3=sqrtS_run3/2.;
    for (int iBin = 1; iBin <= H1D_jetPt_run2_MLPaperFile->GetNbinsX(); iBin++) {
      H1D_jetPt_unfolded_run2Comp_xT[0]->SetBinContent(iBin, dpt_dxt_run2*H1D_jetPt_run2_MLPaperFile->GetBinContent(iBin));
      H1D_jetPt_unfolded_run2Comp_xT[0]->SetBinError(iBin, dpt_dxt_run2*H1D_jetPt_run2_MLPaperFile->GetBinError(iBin));
    }
    for (int iBin = 1; iBin <= H1D_jetPt_unfolded->GetNbinsX(); iBin++) {
      H1D_jetPt_unfolded_run2Comp_xT[1]->SetBinContent(iBin, dpt_dxt_run3*H1D_jetPt_unfolded->GetBinContent(iBin));
      H1D_jetPt_unfolded_run2Comp_xT[1]->SetBinError(iBin, dpt_dxt_run3*H1D_jetPt_unfolded->GetBinError(iBin));
    }

    double fitPtRange_run2[2] = {xtBinningVectorRun2.front(), xtBinningVectorRun2.back()}; //-1 because Tsallis shape only accurate until 120GeV or so
    double fitPtRange_run3[2] = {xtBinningVectorRun3.front(), xtBinningVectorRun3.back()}; //-1 because Tsallis shape only accurate until 120GeV or so
    std::tuple<TF1*, TMatrixDSym, TFitResultPtr> exponentialFitWithLogTransfoResult_xT_run2 = ExponentialFitWithLogTransfo(H1D_jetPt_unfolded_run2Comp_xT[0], fitPtRange_run2);
    std::tuple<TF1*, TMatrixDSym, TFitResultPtr> exponentialFitWithLogTransfoResult_xT_run3 = ExponentialFitWithLogTransfo(H1D_jetPt_unfolded_run2Comp_xT[1], fitPtRange_run3);

    int nBins = 200;
    std::vector<double> xtFineBinningVector = MakeConstantSizeBinning(fitPtRange_run3[0], fitPtRange_run3[1], nBins);
    double* xtFineBinning = &xtFineBinningVector[0];
    TString histName_xT_run2 = "rebinWithFit_xTcomp_run2"+partialUniqueSpecifier;
    std::pair<TH1D*, TGraphErrors*> pairResult_xT_run2 = RebinWithFit(H1D_jetPt_unfolded_run2Comp_xT[0], nBins, xtFineBinning, fitPtRange_run3, histName_xT_run2, exponentialFitWithLogTransfoResult_xT_run2);
    TString histName_xT_run3 = "rebinWithFit_xTcomp_run3"+partialUniqueSpecifier;
    std::pair<TH1D*, TGraphErrors*> pairResult_xT_run3 = RebinWithFit(H1D_jetPt_unfolded_run2Comp_xT[1], nBins, xtFineBinning, fitPtRange_run3, histName_xT_run3, exponentialFitWithLogTransfoResult_xT_run3);

    TGraph_xT_run2Comp_fits.push_back(pairResult_xT_run2.second);
    TGraph_xT_run2Comp_fits.push_back(pairResult_xT_run3.second);

    //ratio
    TH1D* H1D_xT_run2Comp_fits_tgraphConverted[2];
    H1D_xT_run2Comp_fits_tgraphConverted[0] = pairResult_xT_run2.first;
    H1D_xT_run2Comp_fits_tgraphConverted[1] = pairResult_xT_run3.first;
    H1D_xT_ratio_run2Comp_fits = (TH1D*)H1D_xT_run2Comp_fits_tgraphConverted[1]->Clone("H1D_xT_ratio_run2Comp_fits"+partialUniqueSpecifier);
    divideSuccessRun2_xT_fits = H1D_xT_ratio_run2Comp_fits->Divide(H1D_xT_run2Comp_fits_tgraphConverted[0]); // run3/run2

    // ////////// using fit: commented for now as the fits aren't great //////////
    // // make xT binning
    // double binWidth = 0.001;
    // float maxPt = ptWindowDisplay[1];
    // float maxXt = 2*maxPt/sqrtS_run2; // sqrtS_run3 is larger than sqrtS_run2 so maxXtRun3 is smaller than maxXtRun2
    // for (int iBin = 1; iBin <= maxXt/binWidth; iBin++) { //pt bins of 1GeV are travelled
    //   xtBinningVector.push_back(binWidth * iBin);
    // }
    // int nBinsXt = xtBinningVector.size() - 1;
    // double* xtBinning = &xtBinningVector[0];

    // // get fit functions
    // std::tuple<TF1*, TMatrixDSym, TFitResultPtr> tupleFitResult_run2 = TsallisFit(H1D_jetPt_run2_MLPaperFile, nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius], fitPtRange);
    // TF1* TF1_jetPt_run2_fit = std::get<0>(tupleFitResult_run2);
    // TMatrixDSym covMatrix_run2_fit = std::get<1>(tupleFitResult_run2);
    // double fitPtRange_run3[2] = {ptBinsJetsGen[iRadius][0], ptBinsJetsGen[iRadius][nBinPtJetsGen[iRadius]]};
    // std::tuple<TF1*, TMatrixDSym, TFitResultPtr> tupleFitResult_run3 = TsallisFit(H1D_jetPt_unfolded, nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius], fitPtRange);
    // TF1* TF1_jetPt_run3_fit = std::get<0>(tupleFitResult_run3);
    // TMatrixDSym covMatrix_run3_fit = std::get<1>(tupleFitResult_run3);


    // // TF1_jetPt_run2_fit->SetParameters()

    // double parfitFunctionRun2[2];
    // double parfitFunctionRun3[2];
    // TF1_jetPt_run2_fit->GetParameters(&parfitFunctionRun2[0]);
    // TF1_jetPt_run3_fit->GetParameters(&parfitFunctionRun3[0]);

    // // TF1* initial function = new TF1("dNdptFunction_run2", "x*(1+1/([0]*[1])*x)**(-[0])", xtBinning[0], xtBinning[nBinsXt]);
    // // (FoG)'(x)=F'oG(x)*G'(x)
    // // [2] is dpt/dxt
    // // [3] is replacing x with xt=2pt/sqrtS
    // // TF1* dNdxtFunction_run2 = new TF1("dNdxtFunction_run2", "[2]*[3]*x*(1+1/([0]*[1])*[3]*x)**(-[0])", xtBinning[0], xtBinning[nBinsXt]); 
    // TF1* dNdxtFunction_run2 = new TF1("dNdxtFunction_run2", "x*(1+1/([0]*[1])*[2]*x)**(-[0])", xtBinning[0], xtBinning[nBinsXt]); //first [2]*[3] term, with [2] being dpt/dxt=sqrtS/2 and [3] being 2./sqrtS, cancel each other 
    // dNdxtFunction_run2->SetParameters(parfitFunctionRun2[0], parfitFunctionRun2[1], 2./sqrtS_run2); //transfor of pt -> xt in variable: xt=2*pt/sqrtS ; pt=sqrtS/2*xt ; dN/dxt(xt) = dN/dpt*dpt/dxt = dN/dpt(2*pt/sqrtS)*sqrtS/2
    // // TF1* dNdxtFunction_run3 = new TF1("dNdxtFunction_run3", "[2]*[3]*x*(1+1/([0]*[1])*[3]*x)**(-[0])", xtBinning[0], xtBinning[nBinsXt]);
    // TF1* dNdxtFunction_run3 = new TF1("dNdxtFunction_run3", "x*(1+1/([0]*[1])*[2]*x)**(-[0])", xtBinning[0], xtBinning[nBinsXt]);
    // dNdxtFunction_run3->SetParameters(parfitFunctionRun3[0], parfitFunctionRun3[1], 2./sqrtS_run2); //transfor of pt -> xt in variable: xt=2*pt/sqrtS ; pt=sqrtS/2*xt ; dN/dxt(xt) = dN/dpt*dpt/dxt = dN/dpt(2*pt/sqrtS)*sqrtS/2
    // // 2./sqrtS_run2


    // int nRowsNew = 3; // for parameters [0], [1], [2]
    // TMatrixDSym newCovMatrix_run2 = TMatrixDSym(nRowsNew);
    // TMatrixDSym newCovMatrix_run3 = TMatrixDSym(nRowsNew);
    // std::vector<double> initialisationVector;
    // for (int i = 0; i < nRowsNew*nRowsNew; i++) {
    //   initialisationVector.push_back(0);
    // }
    // newCovMatrix_run2.SetMatrixArray(&initialisationVector[0]);
    // newCovMatrix_run2.SetSub(0, 0, covMatrix_run2_fit); // inserts covMatrix_run2_fit as submatrix at row0 column0; other errors are 0 given [2] and[3] are constants
    // newCovMatrix_run3.SetMatrixArray(&initialisationVector[0]);
    // newCovMatrix_run3.SetSub(0, 0, covMatrix_run3_fit); // inserts covMatrix_run2_fit as submatrix at row0 column0; other errors are 0 given [2] and[3] are constants

    // Double_t *pData = newCovMatrix_run2.GetMatrixArray();
    // for (int i = 0; i < nRowsNew*nRowsNew; i++) {
    //   cout << "test i = " << i << ", newCovMatrix_run2[i]" << pData[i] << endl;
    // }

    // double xtRange[2] = {xtBinningVector.front(), xtBinningVector.back()};
    // TGraphErrors* fitFunctionTGraphErrors_run2 = getFunctionTGraphErrorsFromCovMatrix(xtRange, dNdxtFunction_run2, &newCovMatrix_run2, nBinsXt);    
    // TGraphErrors* fitFunctionTGraphErrors_run3 = getFunctionTGraphErrorsFromCovMatrix(xtRange, dNdxtFunction_run3, &newCovMatrix_run3, nBinsXt);    

    // // define xT histograms
    // H1D_jetPt_unfolded_run2Comp_xT[0] = new TH1D("H1D_jetPt_unfolded_run2Comp_xT_run2", "H1D_jetPt_unfolded_run2Comp_xT_run2", nBinsXt, xtBinning);
    // H1D_jetPt_unfolded_run2Comp_xT[1] = new TH1D("H1D_jetPt_unfolded_run2Comp_xT_run3", "H1D_jetPt_unfolded_run2Comp_xT_run3", nBinsXt, xtBinning);
    // H1D_jetPt_unfolded_run2Comp_xT[0]->Sumw2();
    // H1D_jetPt_unfolded_run2Comp_xT[1]->Sumw2();


    // // fill the histograms
    // double H1D_jetPt_unfolded_run2Comp_xT_errorRun2, H1D_jetPt_unfolded_run2Comp_xT_errorRun3;
    // double xtAtCenterOfBin;
    // for (int iBin = 1; iBin <= nBinsXt; iBin++) {
    //   xtAtCenterOfBin = (xtBinning[iBin-1]+xtBinning[iBin])/2;
    //   H1D_jetPt_unfolded_run2Comp_xT[0]->SetBinContent(iBin, dNdxtFunction_run2->Eval(xtAtCenterOfBin));
    //   H1D_jetPt_unfolded_run2Comp_xT[0]->SetBinError(iBin, fitFunctionTGraphErrors_run2->GetErrorY(iBin-1));
    //   H1D_jetPt_unfolded_run2Comp_xT[1]->SetBinContent(iBin, dNdxtFunction_run3->Eval(xtAtCenterOfBin));
    //   H1D_jetPt_unfolded_run2Comp_xT[1]->SetBinError(iBin, fitFunctionTGraphErrors_run3->GetErrorY(iBin-1));
    //   cout << "test fitFunctionTGraphErrors_run2->GetErrorY(iBin) = " << fitFunctionTGraphErrors_run2->GetErrorY(iBin-1) << endl;
    //   cout << "errors of xt graph look way too small (1E-8 or 1E-10)" << endl;
    // }

    // //ratio
    // H1D_jetPt_ratio_run2Comp_xT = (TH1D*)H1D_jetPt_unfolded_run2Comp_xT[1]->Clone("H1D_jetPt_ratio_run2Comp_xT"+partialUniqueSpecifier);
    // divideSuccessRun2_xt = H1D_jetPt_ratio_run2Comp_xT->Divide(H1D_jetPt_unfolded_run2Comp_xT[0]); // run3/run2




    // comparison with run2 results with fits only
    H1D_jetPt_unfolded_run2Comp_fits[0] = (TH1D*)H1D_jetPt_run2_MLPaperFile->Clone("H1D_jetPt_unfolded_run2Comp_fits_run2"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_run2Comp_fits[1] = (TH1D*)H1D_jetPt_unfolded->Clone("H1D_jetPt_unfolded_run2Comp_fits_run3"+partialUniqueSpecifier);

    std::tuple<TF1*, TMatrixDSym, TFitResultPtr> exponentialFitWithLogTransfoResult_run3 = ExponentialFitWithLogTransfo(H1D_jetPt_unfolded, fitPtRange);
    TString histName_run3_fine = "rebinWithFit_run3_fine"+partialUniqueSpecifier;
    TString histName_run2_fine = "rebinWithFit_run3_fine"+partialUniqueSpecifier;
    std::pair<TH1D*, TGraphErrors*> pairResult_run3_fine = RebinWithFit(H1D_jetPt_unfolded, nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], fitPtRange, histName_run3_fine, exponentialFitWithLogTransfoResult_run3);
    std::pair<TH1D*, TGraphErrors*> pairResult_run2_fine = RebinWithFit(H1D_jetPt_run2_MLPaperFile, nBinPtJetsFine[iRadius], ptBinsJetsFine[iRadius], fitPtRange, histName_run2_fine, exponentialFitWithLogTransfoResult_run2);

    TGraph_jetPt_unfolded_run2Comp_fits.push_back(pairResult_run2_fine.second);
    TGraph_jetPt_unfolded_run2Comp_fits.push_back(pairResult_run3_fine.second);
    

    //ratio
    TH1D* H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[2];
    H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[0] = pairResult_run2_fine.first;
    H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[1] = pairResult_run3_fine.first;
    H1D_jetPt_ratio_run2Comp_fits = (TH1D*)H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[1]->Clone("H1D_jetPt_ratio_run2Comp_fits"+partialUniqueSpecifier);
    divideSuccessRun2_fits = H1D_jetPt_ratio_run2Comp_fits->Divide(H1D_jetPt_unfolded_run2Comp_fits_tgraphConverted[0]); // run3/run2
  }

  cout << "comparison with refolded" << endl; 
  if (!useFineBinningTest) {
    Get_Pt_spectrum_bkgCorrected_recBinning(H1D_jetPt_measured, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_bkgCorrected_fineBinning(H1D_jetPt_measured, iDataset, iRadius, options);
  }
  Get_Pt_spectrum_dataUnfoldedThenRefolded(H1D_jetPt_unfoldedThenRefolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options);
  H1D_jetPt_unfolded_refoldedComp[0] = (TH1D*)H1D_jetPt_unfoldedThenRefolded->Clone("H1D_jetPt_refolded_refoldedComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_refoldedComp[1] = (TH1D*)H1D_jetPt_measured->Clone("H1D_jetPt_measured_refoldedComp"+partialUniqueSpecifier);
  H1D_jetPt_ratio_measuredRefolded[0] = (TH1D*)H1D_jetPt_unfoldedThenRefolded->Clone("H1D_jetPt_ratio_refoldedComp"+partialUniqueSpecifier);
  // divideSuccessMeasuredRefolded[0] = H1D_jetPt_ratio_measuredRefolded[0]->Divide(H1D_jetPt_measured);
  divideSuccessMeasuredRefolded = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_measuredRefolded[0], H1D_jetPt_measured);


  if (doClosure_splitMC_mcdUnfoldedVsGen) {
    std::string detRespOption = useFactorisedMatrixInMcdUnfoldedClosure ? ", useIdentityForFluctResp" : "";
    unfoldParameter_mcSplitClosure = Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded_inputSplitClosure, measuredInput_mcSplitInput, iDataset, iRadius, unfoldParameterInput, options+(std::string)detRespOption, splitTestControlMC).first;

    H1D_jetPt_unfolded_mcdSplitClosure[0] = (TH1D*)H1D_jetPt_mcp_mcSplitInput->Clone("H1D_jetPt_mcp_mcSplitClosure"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_mcdSplitClosure[1] = (TH1D*)H1D_jetPt_unfolded_inputSplitClosure->Clone("H1D_jetPt_unfolded_inputSplitClosure"+partialUniqueSpecifier);
    H1D_jetPt_ratio_mcdSplitClosure = (TH1D*)H1D_jetPt_unfolded_inputSplitClosure->Clone("H1D_jetPt_ratio_mcSplitClosure"+partialUniqueSpecifier);
    // divideSuccessMcdSplitClosure = H1D_jetPt_ratio_mcdSplitClosure->Divide(H1D_jetPt_mcp_mcSplitInput);
    divideSuccessMcdSplitClosure = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_mcdSplitClosure, H1D_jetPt_mcp_mcSplitInput);

  }

  if (doClosure_splitMC_mcpFoldedWithFluct) {
    cout << "comparison mcp folded with fluctuations vs mcp" << endl; 
    Get_Pt_spectrum_mcpFoldedWithFluctuations(H1D_jetPt_mcpFolded, iDataset, iRadius, options); // 
    Get_Pt_spectrum_mcp_recBinning(H1D_jetPt_mcp_mcSplitInput_recBinning, iDataset, iRadius, options);

    H1D_jetPt_unfolded_mcpFoldedComp[0] = (TH1D*)H1D_jetPt_mcpFolded->Clone("H1D_jetPt_mcpFoldedComp_mcpFolded"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_mcpFoldedComp[1] = (TH1D*)H1D_jetPt_mcp_mcSplitInput_recBinning->Clone("H1D_jetPt_mcpFoldedComp_mcp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_mcpFoldedMcp = (TH1D*)H1D_jetPt_mcpFolded->Clone("H1D_jetPt_ratio_mcpFoldedMcp"+partialUniqueSpecifier);
    // divideSuccessMcpFoldedMcp = H1D_jetPt_ratio_mcpFoldedMcp->Divide(H1D_jetPt_mcp_mcSplitInput_recBinning);
    divideSuccessMcpFoldedMcp = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_mcpFoldedMcp,H1D_jetPt_mcp_mcSplitInput_recBinning);

    
    // cout << "Integral mcp folded: " << H1D_jetPt_mcpFolded->Integral(1, H1D_jetPt_mcpFolded->GetNbinsX()) << endl;
    // cout << "Integral mcp       : " << H1D_jetPt_mcp_mcSplitInput->Integral(1, H1D_jetPt_mcp_mcSplitInput->GetNbinsX()) << endl;
  
    cout << "comparison mcp folded with fluctuations then unfolded vs mcp" << endl; 
    if (!normGenAndMeasByNEvtsForUnfoldingInput) {
      Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcpFolded2, iDataset, iRadius, options);
    } else{
      Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEnd(H1D_jetPt_mcpFolded2, iDataset, iRadius, options);
    }  
    Get_Pt_spectrum_unfolded(H1D_jetPt_mcpFoldedThenUnfolded, H1D_jetPt_mcpFolded2, iDataset, iRadius, unfoldParameterInput, options+", noKineEff, noPurity, noEff, inputIsMC, inputIsMCPFoldedTest, useIdentityForDetResp"); // input is mcp with fluctuations smearing: there are no fake jets, and the ptbinrange is the gen one so no kine efficiency
    H1D_jetPt_unfolded_mcpFoldedUnfoldedComp[0] = (TH1D*)H1D_jetPt_mcpFoldedThenUnfolded->Clone("H1D_jetPt_mcpFoldedUnfoldedComp_mcpFoldedUnfolded"+partialUniqueSpecifier);
    H1D_jetPt_unfolded_mcpFoldedUnfoldedComp[1] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_mcpFoldedUnfoldedComp_mcp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_mcpFoldedUnfoldedMcp = (TH1D*)H1D_jetPt_mcpFoldedThenUnfolded->Clone("H1D_jetPt_ratio_mcpFoldedUnfoldedMcp"+partialUniqueSpecifier);
    // divideSuccessMcpFoldedUnfoldedMcp = H1D_jetPt_ratio_mcpFoldedUnfoldedMcp->Divide(H1D_jetPt_mcp); // divided by mcp Get_Pt_spectrum_mcp_genBinning
    divideSuccessMcpFoldedUnfoldedMcp = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_mcpFoldedUnfoldedMcp, H1D_jetPt_mcp); // divided by mcp Get_Pt_spectrum_mcp_genBinning
  }

  TString unfoldingCode;
  if (useManualRespMatrixSettingMethod){
    unfoldingCode = "myUnfold";
  } else {
    unfoldingCode = "joUnfold";
  }
  TString unfoldingInfo = (TString)unfoldingMethod+"-k="+Form("%i", unfoldParameter)+"-"+(TString)unfoldingPrior+"-"+unfoldingCode+"-matrixTransfo"+matrixTransformationOrder;

  std::error_code errPDF, errPNG, errEPS;
  CreateDirectoryRecursive((std::string)"pdfFolder/IterationsDump", errPDF);
  CreateDirectoryRecursive((std::string)"pngFolder/IterationsDump", errPNG);
  CreateDirectoryRecursive((std::string)"epsFolder/IterationsDump", errEPS);
  // struct stat st1{};
  // if (stat("pdfFolder/IterationsDump", &st1) == -1) {
  //     mkdir("pdfFolder/IterationsDump", 0700);
  // }
  // struct stat st2{};
  // if (stat("pngFolder/IterationsDump", &st2) == -1) {
  //     mkdir("pngFolder/IterationsDump", 0700);
  // }
  // struct stat st3{};
  // if (stat("epsFolder/ResponseMatrices", &st3) == -1) {
  //     mkdir("epsFolder/ResponseMatrices", 0700);
  // }

  TString textContext = contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), "");

  TString dummyLegend[1] = {""};

  TString* yAxisLabel = texJet_dNdeta;
  TString* yAxisLabelXt = texJet_dNdeta;
  if (normaliseDistribsInComparisonPlots || normaliseUnfoldingResultsAtEnd) { //should probably check if having both on doesn't lead to double normalisation
    yAxisLabel = texJet_d2Ndptdeta_EventNorm;
    yAxisLabelXt = texJet_d2Ndxtdeta_EventNorm;
  }

  TString pdfTitleBase = (TString)"IterationsDump/jet_"+unfoldingInfo;//+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_unfolded_";
  // std::array<std::array<float, 2>, 2> drawnWindow = {{{ptWindowDisplay[0], ptWindowDisplay[1]}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

  // comparison with measured
  TString unfoldedMeasuredCompLegend[2] = {"measured (gen binning)", "unfolded data"};
  TString* pdfName_measuredComp = new TString(pdfTitleBase+"_measuredComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_measuredComp, unfoldedMeasuredCompLegend, 2, textContext, pdfName_measuredComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (divideSuccessMeasured){
    TString* pdfName_ratio_measured = new TString(pdfTitleBase+"_measuredComp_ratio");
    Draw_TH1_Histogram(H1D_jetPt_ratio_measured, textContext, pdfName_ratio_measured, texPtX, texRatioUnfoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");
    // TString* pdfName_ratio_measured_zoom = new TString(pdfTitleBase+"_measuredComp_ratio_zoom");
    // Draw_TH1_Histogram(H1D_jetPt_ratio_measured, textContext, pdfName_ratio_measured_zoom, texPtX, texRatioUnfoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine,zoomToOneMedium2");
  }

  // comparison with mcp truth
  TString unfoldedTruthCompLegend[2] = {"mcp truth", "unfolded data"};
  TString* pdfName_mcpComp = new TString(pdfTitleBase+"_mcpComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpComp, unfoldedTruthCompLegend, 2, textContext, pdfName_mcpComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (divideSuccessMcp){
    TString* pdfName_ratio_mcp = new TString(pdfTitleBase+"_mcpComp_ratio");
    Draw_TH1_Histogram(H1D_jetPt_ratio_mcp, textContext, pdfName_ratio_mcp, texPtX, texRatioUnfoldedMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    TString* pdfName_ratio_mcp_zoom = new TString(pdfTitleBase+"_mcpComp_ratio_zoom");
    Draw_TH1_Histogram(H1D_jetPt_ratio_mcp, textContext, pdfName_ratio_mcp_zoom, texPtX, texRatioUnfoldedMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine,zoomToOneMedium2");
  }


  // comparison with refolded
  TString unfoldedRefoldedCompLegend[2] = {"refolded manually", "measured"};
  TString* pdfName_refoldedComp = new TString(pdfTitleBase+"_RefoldedComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_refoldedComp, unfoldedRefoldedCompLegend, 2, textContext, pdfName_refoldedComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (divideSuccessMeasuredRefolded) {
    TString* pdfName_ratio_refoldedComp = new TString(pdfTitleBase+"_RefoldedComp_ratio");
    Draw_TH1_Histogram(H1D_jetPt_ratio_measuredRefolded[0], textContext, pdfName_ratio_refoldedComp, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    TString* pdfName_ratio_refoldedComp_zoom = new TString(pdfTitleBase+"_RefoldedComp_ratio_zoom");
    Draw_TH1_Histogram(H1D_jetPt_ratio_measuredRefolded[0], textContext, pdfName_ratio_refoldedComp_zoom, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine,zoomToOneMedium2");
  }


  // comparison with Run 2
  if (comparePbPbWithRun2) {
    TString unfoldedRun2CompLegend_fitRebin[3] = {"unfolded Run3", "unfolded Run2 ML rebinned", "unfolded Run2 ML initial"};
    TString* pdfName_run2Comp_fitRebin = new TString(pdfTitleBase+"_run2Comp_fitRebin");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_run2Comp_fitRebin, unfoldedRun2CompLegend_fitRebin, 3, textContext, pdfName_run2Comp_fitRebin, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy,fitSingle", TGraph_jetPt_run2_MLPaperFile_fit);
    if (divideSuccessRun2_fitRebin[0] || divideSuccessRun2_fitRebin[1]) {
      TString* pdfName_ratio_run2_fitRebin = new TString(pdfTitleBase+"_run2Comp_fitRebin_ratio");
      Draw_TH1_Histogram(H1D_jetPt_ratio_run2_fitRebin[0], textContext, pdfName_ratio_run2_fitRebin, texPtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge, ratioLine");
    }

    TString unfoldedRun2CompLegend[3] = {"unfolded Run3", "unfolded Run3 rebinned", "unfolded Run2"};
    TString* pdfName_run2Comp = new TString(pdfTitleBase+"_run2Comp");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_run2Comp, unfoldedRun2CompLegend, 3, textContext, pdfName_run2Comp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
    if (divideSuccessRun2[0] || divideSuccessRun2[1]) {
      TString* pdfName_ratio_run2 = new TString(pdfTitleBase+"_run2Comp_ratio");
      Draw_TH1_Histogram(H1D_jetPt_ratio_run2[0], textContext, pdfName_ratio_run2, texPtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge, ratioLine");
    }

    TString unfoldedRun2CompLegend_shapeComp[2] = {"unfolded Run3 rebinned", "unfolded Run2 scaled up to Run 3"};
    TString* pdfName_run2Comp_shapeComp = new TString(pdfTitleBase+"_run2Comp_shapeComp");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_run2Comp_shapeComp, unfoldedRun2CompLegend_shapeComp, 2, textContext, pdfName_run2Comp_shapeComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
    if (divideSuccessRun2[0] || divideSuccessRun2[1]) {
      TString* pdfName_ratio_shapeComp_run2 = new TString(pdfTitleBase+"_run2Comp_shapeComp_ratio");
      Draw_TH1_Histogram(H1D_jetPt_ratio_run2_shapeComp[0], textContext, pdfName_ratio_shapeComp_run2, texPtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge, ratioLine");
    }

    TString unfoldedRun2CompLegend_xtComp[2] = {"Run2", "Run3"};
    TString* pdfName_run2Comp_xtComp = new TString(pdfTitleBase+"_run2Comp_xtComp");
    std::array<std::array<float, 2>, 2> drawnWindowXt = {{{(float)xtBinningVectorRun3.front(), (float)xtBinningVectorRun3.back()}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}
    Draw_TH1_Histograms(H1D_jetPt_unfolded_run2Comp_xT, unfoldedRun2CompLegend_xtComp, 2, textContext, pdfName_run2Comp_xtComp, texXtX, yAxisLabelXt, texCollisionDataInfo, drawnWindowXt, legendPlacementAuto, contextPlacementAuto, "logy,fitCollection", TGraph_xT_run2Comp_fits);
    if (divideSuccessRun2_xT_fits) {
      std::array<std::array<float, 2>, 2> drawnWindowXtRatios = {{{(float)xtBinningVectorRun3.front(), (float)xtBinningVectorRun3.back()}, {0, 2.2}}}; // {{xmin, xmax}, {ymin, ymax}}
      TString* pdfName_ratio_run2_xT_fits = new TString(pdfTitleBase+"_run2Comp_xT_fits_ratio");
      Draw_TH1_Histogram(H1D_xT_ratio_run2Comp_fits, textContext, pdfName_ratio_run2_xT_fits, texXtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindowXtRatios, legendPlacementAuto, contextPlacementAuto, "ratioLine"); //zoomToOneExtraExtra
    }

    TString unfoldedRun2CompLegend_fits[3] = {"unfolded Run2 ML", "unfolded Run3"};
    TString* pdfName_run2Comp_fits = new TString(pdfTitleBase+"_run2Comp_fits");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_run2Comp_fits, unfoldedRun2CompLegend_fits, 2, textContext, pdfName_run2Comp_fits, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy,fitCollection", TGraph_jetPt_unfolded_run2Comp_fits);
    if (divideSuccessRun2_fits) {
      std::array<std::array<float, 2>, 2> drawnWindowFitRatios = {{{ptWindowDisplay[0], ptWindowDisplay[1]}, {0, 2.2}}}; // {{xmin, xmax}, {ymin, ymax}}
      TString* pdfName_ratio_run2_fits = new TString(pdfTitleBase+"_run2Comp_fits_ratio");
      Draw_TH1_Histogram(H1D_jetPt_ratio_run2Comp_fits, textContext, pdfName_ratio_run2_fits, texPtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindowFitRatios, legendPlacementAuto, contextPlacementAuto, "ratioLine");
    }
  }

  if (doClosure_splitMC_mcpFoldedWithFluct) {
    // comparison mcp folded with fluctuations vs mcp
    TString unfoldedMcpFoldedCheckLegend[2] = {"mcp-folded", "mcp"};
    TString* pdfName_McpFoldedCheck = new TString(pdfTitleBase+"_McpFoldedVsMcp");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpFoldedComp, unfoldedMcpFoldedCheckLegend, 2, textContext, pdfName_McpFoldedCheck, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
    if (H1D_jetPt_ratio_mcpFoldedMcp) {
      TString* pdfName_ratio_McpFoldedCheck = new TString(pdfTitleBase+"_McpFoldedVsMcp_ratio");
      Draw_TH1_Histogram(H1D_jetPt_ratio_mcpFoldedMcp, textContext, pdfName_ratio_McpFoldedCheck, texPtX, texRatioMcpFoldedVsMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");
    }

    // comparison mcp folded with fluctuations then unfolded vs mcp
    TString unfoldedMcpFoldedUnfoldedCheckLegend[2] = {"mcp-folded unfolded", "mcp"};
    TString* pdfName_McpFoldedUnfoldedCheck = new TString(pdfTitleBase+"_McpFoldedUnfoldedCheck");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpFoldedUnfoldedComp, unfoldedMcpFoldedUnfoldedCheckLegend, 2, textContext, pdfName_McpFoldedUnfoldedCheck, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
    if (divideSuccessMcpFoldedUnfoldedMcp) {
      TString* pdfName_ratio_McpFoldedUnfoldedCheck = new TString(pdfTitleBase+"_McpFoldedUnfoldedCheck_ratio");
      Draw_TH1_Histogram(H1D_jetPt_ratio_mcpFoldedUnfoldedMcp, textContext, pdfName_ratio_McpFoldedUnfoldedCheck, texPtX, texRatioMcpFoldedUnfoldedMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");
    }
  }

  if (doClosure_splitMC_mcdUnfoldedVsGen) {
    // comparison mcd from controlMC unfolded vs mcp from controlMC
    TString unfoldedMcdClosureCheckLegend[2] = {"mcp", "mcd unfolded"};
    TString* pdfName_McdUnfoldedClosureCheck = new TString(pdfTitleBase+"_McdUnfoldedClosureCheck");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_mcdSplitClosure, unfoldedMcdClosureCheckLegend, 2, textContext, pdfName_McdUnfoldedClosureCheck, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
    if (divideSuccessMcdSplitClosure) {
      TString* pdfName_ratio_McdUnfoldedClosureCheck = new TString(pdfTitleBase+"_McdUnfoldedClosureCheck_ratio");
      Draw_TH1_Histogram(H1D_jetPt_ratio_mcdSplitClosure, textContext, pdfName_ratio_McdUnfoldedClosureCheck, texPtX, texRatioMcdUnfoldedMcd, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");
    }
  }
  // Persist the same normalised spectra used in the plots for numerical inspection.
  for (int bin = 1; bin <= H1D_jetPt_unfolded->GetNbinsX(); ++bin) {
    if (!std::isfinite(H1D_jetPt_unfolded->GetBinContent(bin)) ||
        !std::isfinite(H1D_jetPt_unfolded->GetBinError(bin))) {
      throw std::runtime_error("Non-finite unfolded spectrum");
    }
  }
  // Temporarily disable saving unfolding spectra; uncomment to enable again.
//   TString resultPath = Form("unfolding_%s_R%02d_%s_k%d.root", Datasets[iDataset].Data(),
//                             int(std::lround(10 * arrayRadius[iRadius])), unfoldingMethod, unfoldParameterInput);
//   TFile resultFile(resultPath, "RECREATE");
//   if (resultFile.IsZombie()) throw std::runtime_error("Cannot create unfolding result file");
//   H1D_jetPt_unfolded->Write("unfolded");
//   H1D_jetPt_measured->Write("measured");
//   H1D_jetPt_unfoldedThenRefolded->Write("refolded");
//   H1D_jetPt_ratio_measuredRefolded[0]->Write("refolded_over_measured");
//   H1D_jetPt_mcp->Write("mcp");
//   TNamed configuration("analysis_options", options.c_str());
//   configuration.Write();
//   resultFile.Close();
//   cout << "Saved unfolding spectra: " << resultPath << endl;

}

void Draw_Pt_spectrum_unfolded_parameterVariation_singleDataset(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options) {

  const int nUnfoldIteration = std::floor((unfoldIterationMax - unfoldIterationMin)/step) + 1;

  TH1D* H1D_jetPt_measured;
  TH1D* H1D_jetPt_measured_genBinning;
  TH1D* H1D_jetPt_unfolded[nUnfoldIteration];
  TH1D* H1D_jetPt_unfoldedThenRefolded[nUnfoldIteration];
  TH1D* H1D_jetPt_unfolded_mcpComp[nUnfoldIteration+1];
  TH1D* H1D_jetPt_unfolded_measuredComp[nUnfoldIteration+1];
  TH1D* H1D_jetPt_unfolded_refoldedComp[nUnfoldIteration+1];
  TH1D* H1D_jetPt_mcp;
  TH1D* H1D_jetPt_ratio_mcp[nUnfoldIteration];
  TH1D* H1D_jetPt_ratio_measured[nUnfoldIteration];
  TH1D* H1D_jetPt_ratio_measuredRefolded[nUnfoldIteration];

  bool divideSuccessMcp[nUnfoldIteration];
  bool divideSuccessMeasured[nUnfoldIteration];
  bool divideSuccessMeasuredRefolded[nUnfoldIteration];
  TString partialUniqueSpecifier;


  if (!useFineBinningTest) {
    Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp, iDataset, iRadius, options);
  } else {
    Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp, iDataset, iRadius, options);
  }
  // H1D_jetPt_mcp = (TH1D*)H1D_jetPt_mcp_defaultBin->Rebin(nBinPtJetsGen[iRadius],"jetPt_mcp_rebinned_"+RadiusLegend[iRadius]+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsJetsGen[iRadius]);

  partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);
  Get_Pt_spectrum_bkgCorrected_recBinning(H1D_jetPt_measured, iDataset, iRadius, options);
  Get_Pt_spectrum_bkgCorrected_genBinning(H1D_jetPt_measured_genBinning, iDataset, iRadius, options);

  int unfoldParameterInput = 0;

  TH1D* measuredInput;
  if (!normGenAndMeasByNEvtsForUnfoldingInput) {
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options); 
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options);
    }
  } else{
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
    }
  }

  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){

    // unfoldParameterInput = iUnfoldIteration * step + unfoldIterationMin; 
    unfoldParameterInput = unfoldIterationMax - iUnfoldIteration * step; 
    // unfoldingIterationLegend[iUnfoldIteration] = unfoldParameterInput;

    cout << "((((((((((((()))))))))))))" << endl;
    cout << "Iteration "<< iUnfoldIteration << endl;
    cout << "((((((((((((()))))))))))))" << endl;
    Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iUnfoldIteration], measuredInput, iDataset, iRadius, unfoldParameterInput, options);

    // comparison with measured
    H1D_jetPt_unfolded_measuredComp[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfolded[iUnfoldIteration]->Clone("H1D_jetPt_unfolded_measuredComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_measured[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfolded[iUnfoldIteration]->Clone("H1D_jetPt_ratio_measured"+partialUniqueSpecifier);
    // divideSuccessMeasured[iUnfoldIteration] = H1D_jetPt_ratio_measured[iUnfoldIteration]->Divide(H1D_jetPt_measured_genBinning);
    divideSuccessMeasured[iUnfoldIteration] = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_measured[iUnfoldIteration], H1D_jetPt_measured_genBinning); 

    // comparison with mcp truth
    H1D_jetPt_unfolded_mcpComp[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfolded[iUnfoldIteration]->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_mcp[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfolded[iUnfoldIteration]->Clone("H1D_jetPt_ratio_mcp"+partialUniqueSpecifier);
    divideSuccessMcp[iUnfoldIteration] = H1D_jetPt_ratio_mcp[iUnfoldIteration]->Divide(H1D_jetPt_mcp);


    // comparison with refolded
    Get_Pt_spectrum_dataUnfoldedThenRefolded(H1D_jetPt_unfoldedThenRefolded[iUnfoldIteration], measuredInput, iDataset, iRadius, unfoldParameterInput, options);
    H1D_jetPt_unfolded_refoldedComp[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfoldedThenRefolded[iUnfoldIteration]->Clone("H1D_jetPt_refolded_refoldedComp"+partialUniqueSpecifier);
    H1D_jetPt_ratio_measuredRefolded[iUnfoldIteration] = (TH1D*)H1D_jetPt_unfoldedThenRefolded[iUnfoldIteration]->Clone("H1D_jetPt_ratio_refoldedComp"+partialUniqueSpecifier);
    // divideSuccessMeasuredRefolded[iUnfoldIteration] = H1D_jetPt_ratio_measuredRefolded[iUnfoldIteration]->Divide(H1D_jetPt_measured);
    divideSuccessMeasuredRefolded[iUnfoldIteration] = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_measuredRefolded[iUnfoldIteration], H1D_jetPt_measured);

  }
  H1D_jetPt_unfolded_measuredComp[nUnfoldIteration] = (TH1D*)H1D_jetPt_measured_genBinning->Clone("H1D_jetPt_measured_genBinning_measuredComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_mcpComp[nUnfoldIteration] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_unfolded_mcpComp"+partialUniqueSpecifier);
  H1D_jetPt_unfolded_refoldedComp[nUnfoldIteration] = (TH1D*)H1D_jetPt_measured->Clone("H1D_jetPt_measured_refoldedComp"+partialUniqueSpecifier);


  TString unfoldingInfo = (TString)unfoldingMethod+"_"+(TString)unfoldingPrior+"_kmax="+Form("%i", unfoldIterationMax);

  TString unfoldingIterationLegend[nUnfoldIteration+1]; IterationLegend(unfoldingIterationLegend, unfoldIterationMin, unfoldIterationMax, step);

  
  partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  TString textContext = contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), "");

  TString* pdfName = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo);

  TString* yAxisLabel = texCount;
  if (normaliseDistribsInComparisonPlots || normaliseUnfoldingResultsAtEnd) { //should probably check if having both on doesn't lead to double normalisation
    yAxisLabel = texJet_d2Ndptdeta_EventNorm;
  }
  Draw_TH1_Histograms(H1D_jetPt_unfolded, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");


    // comparison with measured
  // TString unfoldedMeasuredCompLegend[2] = {"unfolded data", "measured (gen binning)"};
  unfoldingIterationLegend[nUnfoldIteration] = (TString)"measured";
  TString* pdfName_measuredComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_measuredComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_measuredComp, unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_measuredComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");

  bool divideSuccessMeasured_boolsum = true;
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    if (!(divideSuccessMeasured[iUnfoldIteration])) {
      divideSuccessMeasured_boolsum = false;
    }
  }
  if (divideSuccessMeasured_boolsum){
    TString* pdfName_ratio_measured = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioMeasured");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measured, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_measured, texPtX, texRatioUnfoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");
  }


    // comparison with mcp truth
  // TString unfoldedTruthCompLegend[2] = {"unfolded data", "mcp truth"};
  unfoldingIterationLegend[nUnfoldIteration] = (TString)"mcp";
  TString* pdfName_mcpComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_mcpComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpComp, unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_mcpComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");

  bool divideSuccessMcp_boolsum = true;
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    if (!(divideSuccessMcp[iUnfoldIteration])) {
      divideSuccessMcp_boolsum = false;
    }
  }
  if (divideSuccessMcp_boolsum){
    TString* pdfName_ratio_mcp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioMcp");
    Draw_TH1_Histograms(H1D_jetPt_ratio_mcp, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_mcp, texPtX, texRatioUnfoldedMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
  }

  // comparison with refolded
  // TString unfoldedRefoldedCompLegend[2] = {"refolded", "measured"};
  unfoldingIterationLegend[nUnfoldIteration] = (TString)"measured";
  for(int iIteration = 0; iIteration < nUnfoldIteration; iIteration++){
    unfoldingIterationLegend[iIteration] += (TString)" refolded";
  }
  TString* pdfName_refoldedComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_RefoldedComp");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_refoldedComp, unfoldingIterationLegend, nUnfoldIteration+1, textContext, pdfName_refoldedComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");

  bool divideSuccessRefoldedComp_boolsum = true;
  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    if (!(divideSuccessMeasuredRefolded[iUnfoldIteration])) {
      divideSuccessRefoldedComp_boolsum = false;
    }
  }
  if (divideSuccessRefoldedComp_boolsum){
    TString* pdfName_ratio_refoldedComp = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioRefoldedUnfolded");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measuredRefolded, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_refoldedComp, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    TString* pdfName_ratio_refoldedComp_zoom = new TString("jet_"+jetType[iJetType]+"_"+jetLevel[iJetLevel]+"_"+partialUniqueSpecifier+"_Pt_unfolded_"+unfoldingInfo+"_ratioRefoldedUnfolded_zoom");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measuredRefolded, unfoldingIterationLegend, nUnfoldIteration, textContext, pdfName_ratio_refoldedComp_zoom, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine,zoomToOneMedium2");
  }
}

void Draw_Pt_spectrum_unfolded_datasetComparison(int iRadius, int unfoldParameterInput, std::string options) {
  bool splitTestControlMC = true;

  TH1D* H1D_jetPt_measured[nDatasets];
  TH1D* H1D_jetPt_measured_genBinning[nDatasets];
  TH1D* H1D_jetPt_unfolded[nDatasets];
  TH1D* H1D_jetPt_unfolded_ratio_datasets[nDatasets];
  TH1D* H1D_jetPt_unfoldedThenRefolded[nDatasets];
  // TH1D* H1D_jetPt_mcpFolded[nDatasets];
  // TH1D* H1D_jetPt_mcpFolded2[nDatasets];
  // TH1D* H1D_jetPt_mcpFoldedThenUnfolded[nDatasets];
  // TH1D* H1D_jetPt_unfolded_mcpComp[2][nDatasets];
  // TH1D* H1D_jetPt_unfolded_measuredComp[2][nDatasets];
  TH1D* H1D_jetPt_unfolded_refoldedComp[nDatasets];
  // TH1D* H1D_jetPt_unfolded_mcpFoldedComp[2][nDatasets];
  // TH1D* H1D_jetPt_unfolded_mcpFoldedUnfoldedComp[2][nDatasets];
  TH1D* H1D_jetPt_mcp[nDatasets];
  TH1D* H1D_jetPt_mcp_recBinControl[nDatasets];
  TH1D* H1D_jetPt_ratio_mcp[nDatasets];
  TH1D* H1D_jetPt_ratio_measured[nDatasets];
  TH1D* H1D_jetPt_ratio_measuredRefolded[nDatasets];

  TH1D* H1D_jetPt_unfolded_run2Comp[nDatasets+1];
  TH1D* H1D_jetPt_ratio_run2[nDatasets];
  
  // TH1D* H1D_jetPt_ratio_mcpFoldedMcp[nDatasets];
  // TH1D* H1D_jetPt_ratio_mcpFoldedUnfoldedMcp[nDatasets];


  bool divideSuccessDatasets[nDatasets];
  bool divideSuccessMcp[nDatasets];
  bool divideSuccessMeasured[nDatasets];
  bool divideSuccessMeasuredRefolded[nDatasets];
  bool divideSuccessRun2[nDatasets];
  // bool divideSuccessMcpFoldedMcp[nDatasets];
  // bool divideSuccessMcpFoldedUnfoldedMcp[nDatasets];

  TString partialUniqueSpecifier = (TString)"datasetComparison_R="+Form("%.1f",arrayRadius[iRadius]);
  TString datasetNameSpecifier[nDatasets];

  int unfoldParameter[nDatasets];


  // RUN 2 settings
  TH1D* H1D_jetPt_run2_MLPaperFile = new TH1D("H1D_jetPt_run2_MLPaperFile", "H1D_jetPt_run2_MLPaperFile", nBinPtJetsGen_run2[iRadius], ptBinsJetsGen_run2[iRadius]);
  if (comparePbPbWithRun2) {
    TGraph* Graph_jetPt_run2_MLPaperFile;
    double Ncoll;
    if (centralityRange[0] == 00 && centralityRange[1] == 10) {
      // Ncoll = (1780.9+1387.0)/2; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
      Ncoll = (1956+1722+1521+1346)/4; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
    } else if (centralityRange[0] == 50 && centralityRange[1] == 70) {
      // Ncoll = (103.7+46.1)/2; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
      Ncoll = (89.8+39.8)/2; // https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf in Run 3, https://alice-notes.web.cern.ch/system/files/notes/analysis/453/2017-Sep-26-analysis_note-ALICE_analysis_note.pdf in Run 2
    } else {
      cout << "comparison with run2: Ncoll hasn't been calculated for this centrality interval" << endl;
    }
    double sigmaNN = 67.6; // value for sqrt(s) = 5.02 TeV https://arxiv.org/abs/1710.07098
    double T_AA = Ncoll / sigmaNN;
    Graph_jetPt_run2_MLPaperFile = ((TGraph*)((TDirectoryFile*)file_O2Analysis_run2ComparisonFileMLPaper->Get("Figure 3a top R020"))->FindObjectAny("Graph1D_y1")); // https://doi.org/10.1016/j.physletb.2023.138412
    // H1D_jetPt_run2_MLPaperFile = (TH1D*)((TH1D*)(file_O2Analysis_run2ComparisonFileMLPaper->Get("Figure 3a top R020"))->FindObject("Graph1D_y1"))->Clone("H1D_jetPt_run2_MLPaperFile");
    int Ngraph = Graph_jetPt_run2_MLPaperFile->GetN();
    for (int i=0; i < Ngraph; ++i) // setting bin contents to the TGraph values
    {
      double x,y;
      Graph_jetPt_run2_MLPaperFile->GetPoint(i, x, y);
      H1D_jetPt_run2_MLPaperFile->Fill(x, y); // uncertainties are of course screwed up
      int iHist = H1D_jetPt_run2_MLPaperFile->GetXaxis()->FindBin(x);
      H1D_jetPt_run2_MLPaperFile->SetBinError(iHist, Graph_jetPt_run2_MLPaperFile->GetErrorY(i));
    }
    H1D_jetPt_run2_MLPaperFile->Scale(T_AA);
    // now the spectre from the file is 1/N d2N/dpTdeta, instead of 1/T_AA 1/N d2N/dpTdeta
    H1D_jetPt_unfolded_run2Comp[nDatasets] = (TH1D*)H1D_jetPt_run2_MLPaperFile->Clone("H1D_jetPt_unfolded_run2Comp_run2Rebin_Run2"+partialUniqueSpecifier);
  }

  for (int iDataset = 0; iDataset < nDatasets; ++iDataset) {
    // getting inputs to unfolding
    if (!useFineBinningTest) {
      Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp[iDataset], iDataset, iRadius, options);
      if (doClosure_splitMC_mcdUnfoldedVsGen) {
        Get_Pt_spectrum_mcp_genBinning(H1D_jetPt_mcp_recBinControl[iDataset], iDataset, iRadius, options, splitTestControlMC);
      }
    } else {
      Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp[iDataset], iDataset, iRadius, options);
      if (doClosure_splitMC_mcdUnfoldedVsGen) {
        Get_Pt_spectrum_mcp_fineBinning(H1D_jetPt_mcp_recBinControl[iDataset], iDataset, iRadius, options, splitTestControlMC);
      }
    }

    TH1D* measuredInput[nDatasets];
    if (!normGenAndMeasByNEvtsForUnfoldingInput) {
      Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput[iDataset], iDataset, iRadius, options); 
      if (useFineBinningTest) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput[iDataset], iDataset, iRadius, options);
      }
    } else{
      Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput[iDataset], iDataset, iRadius, options);
      if (useFineBinningTest) {
        Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput[iDataset], iDataset, iRadius, options);
      }
    }

    // doing the unfolding
    unfoldParameter[iDataset] = Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iDataset], measuredInput[iDataset], iDataset, iRadius, unfoldParameterInput, options).first;
    // TH1D* H1D_jetPt_unfolded2 = (TH1D*)H1D_jetPt_unfolded->Clone(H1D_jetPt_unfolded->GetName()+(TString)"H1D_jetPt_unfolded2");
    H1D_jetPt_unfolded_ratio_datasets[iDataset] = (TH1D*)H1D_jetPt_unfolded[iDataset]->Clone("H1D_jetPt_unfolded_ratio_datasets"+partialUniqueSpecifier+datasetNameSpecifier[iDataset]);
    divideSuccessDatasets[iDataset] = H1D_jetPt_unfolded_ratio_datasets[iDataset]->Divide(H1D_jetPt_unfolded[0]);
    // Creating to-be-plotted histograms
    datasetNameSpecifier[iDataset] = "_"+DatasetsNames[iDataset]+Form("%.1d",iDataset);


    cout << "comparison with measured" << endl; 
    if (!useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_genBinning(H1D_jetPt_measured_genBinning[iDataset], iDataset, iRadius, options);
    } else {
      Get_Pt_spectrum_bkgCorrected_fineBinning(H1D_jetPt_measured_genBinning[iDataset], iDataset, iRadius, options);
    }
    H1D_jetPt_ratio_measured[iDataset] = (TH1D*)H1D_jetPt_unfolded[iDataset]->Clone("H1D_jetPt_ratio_measured"+partialUniqueSpecifier+datasetNameSpecifier[iDataset]);
    // divideSuccessMeasured[iDataset] = H1D_jetPt_ratio_measured[iDataset]->Divide(H1D_jetPt_measured_genBinning[iDataset]);
    divideSuccessMeasured[iDataset] = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_measured[iDataset], H1D_jetPt_measured_genBinning[iDataset]);

    cout << "comparison with mcp truth" << endl; 
    H1D_jetPt_ratio_mcp[iDataset] = (TH1D*)H1D_jetPt_unfolded[iDataset]->Clone("H1D_jetPt_ratio_mcp"+partialUniqueSpecifier+datasetNameSpecifier[iDataset]);
    divideSuccessMcp[iDataset] = H1D_jetPt_ratio_mcp[iDataset]->Divide(H1D_jetPt_mcp[iDataset]);

    cout << "comparison with refolded" << endl; 
    if (!useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_recBinning(H1D_jetPt_measured[iDataset], iDataset, iRadius, options);
    } else {
      Get_Pt_spectrum_bkgCorrected_fineBinning(H1D_jetPt_measured[iDataset], iDataset, iRadius, options);
    }
    Get_Pt_spectrum_dataUnfoldedThenRefolded(H1D_jetPt_unfoldedThenRefolded[iDataset], measuredInput[iDataset], iDataset, iRadius, unfoldParameterInput, options);
    H1D_jetPt_unfolded_refoldedComp[iDataset] = (TH1D*)H1D_jetPt_unfoldedThenRefolded[iDataset]->Clone("H1D_jetPt_refolded_refoldedComp"+partialUniqueSpecifier+datasetNameSpecifier[iDataset]);
    H1D_jetPt_ratio_measuredRefolded[iDataset] = (TH1D*)H1D_jetPt_unfoldedThenRefolded[iDataset]->Clone("H1D_jetPt_ratio_refoldedComp"+partialUniqueSpecifier+datasetNameSpecifier[iDataset]);
    // divideSuccessMeasuredRefolded[iDataset] = H1D_jetPt_ratio_measuredRefolded[iDataset]->Divide(H1D_jetPt_measured[iDataset]);
    divideSuccessMeasuredRefolded[iDataset] = DivideWithCorrelatedErrors_simpleMax(H1D_jetPt_ratio_measuredRefolded[iDataset], H1D_jetPt_measured[iDataset]);

    cout << "comparison with run2" << endl; 
    if (comparePbPbWithRun2) {
      // // comparison with run2 results rebinned using a fit (errors look way underestimated; tsallis function not great aboe 100+ GeV; where should one eval the function inside a bin? probably not just the center)
      // H1D_jetPt_unfolded_run2Comp_fitRebin[0][iDataset] = (TH1D*)H1D_jetPt_unfolded[iDataset]->Clone("H1D_jetPt_unfolded_run2Comp_fitRebin"+partialUniqueSpecifier);
      // double fitPtRange[2] = {ptBinsJetsGen_run2[iRadius][0], ptBinsJetsGen_run2[iRadius][nBinPtJetsGen_run2[iRadius]]};
      // std::pair<TH1D*, TF1*> pairResult = RebinWithTsallisFit(H1D_jetPt_run2_MLPaperFile, nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius], fitPtRange);
      // H1D_jetPt_run2_MLPaperFile_rebinned = pairResult.first;
      // TF1_jetPt_run2_MLPaperFile_fit[0] = pairResult.second;
      // H1D_jetPt_unfolded_run2Comp_fitRebin[1][iDataset] = (TH1D*)H1D_jetPt_run2_MLPaperFile_rebinned->Clone("H1D_jetPt_unfolded_run2_rebinned_fitRebin"+partialUniqueSpecifier);
      // H1D_jetPt_unfolded_run2Comp_fitRebin[2][iDataset] = (TH1D*)H1D_jetPt_run2_MLPaperFile->Clone("H1D_jetPt_unfolded_run2_fitRebin"+partialUniqueSpecifier);
      // // H1D_jetPt_unfolded_run2Comp[2] = (TH1D*)H1D_jetPt_run2_HannaBossiLauraFile->Clone("H1D_jetPt_unfolded_run2Comp_HannaBossiLauraFile"+partialUniqueSpecifier);


      // H1D_jetPt_ratio_run2_fitRebin[0][iDataset] = (TH1D*)H1D_jetPt_run2_MLPaperFile_rebinned->Clone("H1D_jetPt_ratio_run2_MLPaperFile_fitRebin"+partialUniqueSpecifier);
      // // H1D_jetPt_ratio_run2[1] = (TH1D*)H1D_jetPt_run2_HannaBossiLauraFile->Clone("H1D_jetPt_ratio_run2_HannaBossiLauraFile"+partialUniqueSpecifier);
      // divideSuccessRun2_fitRebin[0][iDataset] = H1D_jetPt_ratio_run2_fitRebin[0][iDataset]->Divide(H1D_jetPt_unfolded[iDataset]);
      // // divideSuccessRun2[1] = H1D_jetPt_ratio_run2[1]->Divide(H1D_jetPt_unfolded);


      // comparison with run2 results by rebinning Run3 into Run2 bins
      H1D_jetPt_unfolded_run2Comp[iDataset] = (TH1D*)H1D_jetPt_unfolded[iDataset]->Rebin(nBinPtJetsGen_run2[iRadius],"H1D_jetPt_unfolded_run2Comp_run2Rebin"+partialUniqueSpecifier+datasetNameSpecifier[iDataset], ptBinsJetsGen_run2[iRadius]);
      int scalingFactorRebin[9] = {2, 2, 2, 2, 2, 3, 3, 1, 1}; //width of run 2 histogram
      for (auto i = 1; i <= H1D_jetPt_unfolded_run2Comp[iDataset]->GetNbinsX(); i++) {
        H1D_jetPt_unfolded_run2Comp[iDataset]->SetBinContent(i, 1./scalingFactorRebin[i-1]*H1D_jetPt_unfolded_run2Comp[iDataset]->GetBinContent(i));
        H1D_jetPt_unfolded_run2Comp[iDataset]->SetBinError(i, 1./scalingFactorRebin[i-1]*H1D_jetPt_unfolded_run2Comp[iDataset]->GetBinError(i));
      }
      H1D_jetPt_ratio_run2[iDataset] = (TH1D*)H1D_jetPt_unfolded_run2Comp[iDataset]->Clone("H1D_jetPt_ratio_run2_MLPaperFile"+partialUniqueSpecifier+datasetNameSpecifier[iDataset]);
      divideSuccessRun2[iDataset] = H1D_jetPt_ratio_run2[iDataset]->Divide(H1D_jetPt_run2_MLPaperFile);
    }

    // if (doClosure_splitMC_mcpFoldedWithFluct) {
    //   cout << "comparison mcp folded with fluctuations vs mcp" << endl; 
    //   Get_Pt_spectrum_mcpFoldedWithFluctuations(H1D_jetPt_mcpFolded, iDataset, iRadius, options);
    //   H1D_jetPt_unfolded_mcpFoldedComp[0] = (TH1D*)H1D_jetPt_mcpFolded->Clone("H1D_jetPt_mcpFoldedComp_mcpFolded"+partialUniqueSpecifier);
    //   H1D_jetPt_unfolded_mcpFoldedComp[1] = (TH1D*)H1D_jetPt_mcp_recBinControl->Clone("H1D_jetPt_mcpFoldedComp_mcp"+partialUniqueSpecifier);
    //   H1D_jetPt_ratio_mcpFoldedMcp = (TH1D*)H1D_jetPt_mcpFolded->Clone("H1D_jetPt_ratio_mcpFoldedMcp"+partialUniqueSpecifier);
    //   divideSuccessMcpFoldedMcp = H1D_jetPt_ratio_mcpFoldedMcp->Divide(H1D_jetPt_mcp_recBinControl);

    //   // cout << "Integral mcp folded: " << H1D_jetPt_mcpFolded->Integral(1, H1D_jetPt_mcpFolded->GetNbinsX()) << endl;
    //   // cout << "Integral mcp       : " << H1D_jetPt_mcp_recBinControl->Integral(1, H1D_jetPt_mcp_recBinControl->GetNbinsX()) << endl;
    
    //   cout << "comparison mcp folded with fluctuations then unfolded vs mcp" << endl; 
    //   if (!normGenAndMeasByNEvtsForUnfoldingInput) {
    //     Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEndAndEvtNorm(H1D_jetPt_mcpFolded2, iDataset, iRadius, options);
    //   } else{
    //     Get_Pt_spectrum_mcpFoldedWithFluctuations_preWidthScalingAtEnd(H1D_jetPt_mcpFolded2, iDataset, iRadius, options);
    //   }  
    //   Get_Pt_spectrum_unfolded(H1D_jetPt_mcpFoldedThenUnfolded, H1D_jetPt_mcpFolded2, iDataset, iRadius, unfoldParameterInput, options+", noKineEff, noPurity, noEff, inputIsMC, inputIsMCPFoldedTest"); // input is mcp with fluctuations smearing: there are no fake jets, and the ptbinrange is the gen one so no kine efficiency
    //   H1D_jetPt_unfolded_mcpFoldedUnfoldedComp[0] = (TH1D*)H1D_jetPt_mcpFoldedThenUnfolded->Clone("H1D_jetPt_mcpFoldedUnfoldedComp_mcpFoldedUnfolded"+partialUniqueSpecifier);
    //   H1D_jetPt_unfolded_mcpFoldedUnfoldedComp[1] = (TH1D*)H1D_jetPt_mcp->Clone("H1D_jetPt_mcpFoldedUnfoldedComp_mcp"+partialUniqueSpecifier);
    //   H1D_jetPt_ratio_mcpFoldedUnfoldedMcp = (TH1D*)H1D_jetPt_mcpFoldedThenUnfolded->Clone("H1D_jetPt_ratio_mcpFoldedUnfoldedMcp"+partialUniqueSpecifier);
    //   divideSuccessMcpFoldedUnfoldedMcp = H1D_jetPt_ratio_mcpFoldedUnfoldedMcp->Divide(H1D_jetPt_mcp);
    // }
  }

  TString unfoldingCode;
  if (useManualRespMatrixSettingMethod){
    unfoldingCode = "myUnfold";
  } else {
    unfoldingCode = "joUnfold";
  }
  TString unfoldingInfo = (TString)unfoldingMethod+"-k="+Form("%i", unfoldParameterInput)+"-"+(TString)unfoldingPrior+"-"+unfoldingCode+"-matrixTransfo"+matrixTransformationOrder;

  std::error_code errPDF, errPNG, errEPS;
  CreateDirectoryRecursive((std::string)"pdfFolder/IterationsDump", errPDF);
  CreateDirectoryRecursive((std::string)"pngFolder/IterationsDump", errPNG);
  CreateDirectoryRecursive((std::string)"epsFolder/IterationsDump", errEPS);
  // struct stat st1{};
  // if (stat("pdfFolder/IterationsDump", &st1) == -1) {
  //     mkdir("pdfFolder/IterationsDump", 0700);
  // }
  // struct stat st2{};
  // if (stat("pngFolder/IterationsDump", &st2) == -1) {
  //     mkdir("pngFolder/IterationsDump", 0700);
  // }
  // struct stat st3{};
  // if (stat("epsFolder/ResponseMatrices", &st3) == -1) {
  //     mkdir("epsFolder/ResponseMatrices", 0700);
  // }

  TString textContext = contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextJetRadius(arrayRadius[iRadius]), "");

  TString dummyLegend[1] = {""};

  TString* yAxisLabel = texCount;
  if (normaliseDistribsInComparisonPlots || normaliseUnfoldingResultsAtEnd) { //should probably check if having both on doesn't lead to double normalisation
    yAxisLabel = texJet_d2Ndptdeta_EventNorm;
  }

  TString pdfTitleBase = (TString)"IterationsDump/jet_DatasetComp_"+unfoldingInfo;//+Datasets[iDataset]+DatasetsNames[iDataset]+"_R="+Form("%.1f", arrayRadius[iRadius])+"_Pt_unfolded_";
  // std::array<std::array<float, 2>, 2> drawnWindow = {{{ptWindowDisplay[0], ptWindowDisplay[1]}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

  TString* pdfName_unfolded = new TString(pdfTitleBase+"_unfoldedSpectrum");
  Draw_TH1_Histograms(H1D_jetPt_unfolded, DatasetsNames, nDatasets, textContext, pdfName_unfolded, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (std::all_of(std::begin(divideSuccessDatasets), std::end(divideSuccessDatasets), [](bool booleanEntry) {return booleanEntry;})){ // checks all entries of divideSuccessMeasured are true
    TString* pdfName_unfolded_ratio = new TString(pdfTitleBase+"_unfoldedSpectrum_ratio");
    Draw_TH1_Histograms(H1D_jetPt_unfolded_ratio_datasets, DatasetsNames, nDatasets, textContext, pdfName_unfolded_ratio, texPtX, texRatioDatasets, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    // TString* pdfName_ratio_measured_zoom = new TString(pdfTitleBase+"_ratioToMeasured_zoom");
    // Draw_TH1_Histograms(H1D_jetPt_ratio_measured, DatasetsNames, nDatasets, textContext, pdfName_unfolded_ratio, texPtX, texRatioUnfoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine,zoomToOneMedium2");
  }

    // comparison with raw measured
  if (std::all_of(std::begin(divideSuccessMeasured), std::end(divideSuccessMeasured), [](bool booleanEntry) {return booleanEntry;})){ // checks all entries of divideSuccessMeasured are true
    TString* pdfName_ratio_measured = new TString(pdfTitleBase+"_ratioToMeasured");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measured, DatasetsNames, nDatasets, textContext, pdfName_ratio_measured, texPtX, texRatioUnfoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    TString* pdfName_ratio_measured_zoom = new TString(pdfTitleBase+"_ratioToMeasured_zoom");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measured, DatasetsNames, nDatasets, textContext, pdfName_ratio_measured_zoom, texPtX, texRatioUnfoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine,zoomToOneMedium2");
  }

    // comparison with mcp truth
  if (std::all_of(std::begin(divideSuccessMcp),std::end(divideSuccessMcp), [](bool booleanEntry) {return booleanEntry;})){ // checks all entries of divideSuccessMcp are true
    TString* pdfName_ratio_mcp = new TString(pdfTitleBase+"_RatioToMcp");
    Draw_TH1_Histograms(H1D_jetPt_ratio_mcp, DatasetsNames, nDatasets, textContext, pdfName_ratio_mcp, texPtX, texRatioUnfoldedMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    TString* pdfName_ratio_mcp_zoom = new TString(pdfTitleBase+"_RatioToMcp_zoom");
    Draw_TH1_Histograms(H1D_jetPt_ratio_mcp, DatasetsNames, nDatasets, textContext, pdfName_ratio_mcp_zoom, texPtX, texRatioUnfoldedMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine,zoomToOneMedium2");
  }


  // comparison with refolded
  TString* pdfName_refoldedComp = new TString(pdfTitleBase+"_refoldedSpectrum");
  Draw_TH1_Histograms(H1D_jetPt_unfolded_refoldedComp, DatasetsNames, nDatasets, textContext, pdfName_refoldedComp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (std::all_of(std::begin(divideSuccessMeasuredRefolded), std::end(divideSuccessMeasuredRefolded), [](bool booleanEntry) {return booleanEntry;})){ // checks all entries of divideSuccessMeasuredRefolded are true
    TString* pdfName_ratio_refoldedComp = new TString(pdfTitleBase+"_refoldedSpetrum_ratioToMeasured");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measuredRefolded, DatasetsNames, nDatasets, textContext, pdfName_ratio_refoldedComp, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine");
    TString* pdfName_ratio_refoldedComp_zoom = new TString(pdfTitleBase+"_refoldedSpetrum_ratioToMeasured_zoom");
    Draw_TH1_Histograms(H1D_jetPt_ratio_measuredRefolded, DatasetsNames, nDatasets, textContext, pdfName_ratio_refoldedComp_zoom, texPtX, texRatioRefoldedMeasured, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge,ratioLine,zoomToOneMedium2");
  }


  // comparison with run2
  TString* pdfName_run2Comp = new TString(pdfTitleBase+"_run2Comp");
  TString DatasetsNamesAppended[nDatasets+1];
  for (int iDataset = 0; iDataset < nDatasets; ++iDataset) {
    DatasetsNamesAppended[iDataset] = DatasetsNames[iDataset];
  }
  DatasetsNamesAppended[nDatasets] = "Run2 ML";

  Draw_TH1_Histograms(H1D_jetPt_unfolded_run2Comp, DatasetsNamesAppended, nDatasets+1, textContext, pdfName_run2Comp, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  if (std::all_of(std::begin(divideSuccessRun2), std::end(divideSuccessRun2), [](bool booleanEntry) {return booleanEntry;})){ // checks all entries of divideSuccessRun2 are true
    TString* pdfName_ratio_run2 = new TString(pdfTitleBase+"_run2Comp_ratio");
    Draw_TH1_Histograms(H1D_jetPt_ratio_run2, DatasetsNames, nDatasets, textContext, pdfName_ratio_run2, texPtX, texRatioRun2Unfolded, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "zoomToOneLarge, ratioLine");
  }

  // if (doClosure_splitMC_mcpFoldedWithFluct) {
  //   // comparison mcp folded with fluctuations vs mcp
  //   TString unfoldedMcpFoldedCheckLegend[2] = {"mcp-folded", "mcp"};
  //   TString* pdfName_McpFoldedCheck = new TString(pdfTitleBase+"_McpFoldedVsMcp");
  //   Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpFoldedComp, unfoldedMcpFoldedCheckLegend, 2, textContext, pdfName_McpFoldedCheck, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  //   if (H1D_jetPt_ratio_mcpFoldedMcp) {
  //     TString* pdfName_ratio_McpFoldedCheck = new TString(pdfTitleBase+"_McpFoldedVsMcp_ratio");
  //     Draw_TH1_Histogram(H1D_jetPt_ratio_mcpFoldedMcp, textContext, pdfName_ratio_McpFoldedCheck, texPtX, texRatioMcpFoldedVsMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");
  //   }

  //   // comparison mcp folded with fluctuations then unfolded vs mcp
  //   TString unfoldedMcpFoldedUnfoldedCheckLegend[2] = {"mcp-folded unfolded", "mcp"};
  //   TString* pdfName_McpFoldedUnfoldedCheck = new TString(pdfTitleBase+"_McpFoldedUnfoldedCheck");
  //   Draw_TH1_Histograms(H1D_jetPt_unfolded_mcpFoldedUnfoldedComp, unfoldedMcpFoldedUnfoldedCheckLegend, 2, textContext, pdfName_McpFoldedUnfoldedCheck, texPtX, yAxisLabel, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "logy");
  //   if (divideSuccessMcpFoldedUnfoldedMcp) {
  //     TString* pdfName_ratio_McpFoldedUnfoldedCheck = new TString(pdfTitleBase+"_McpFoldedUnfoldedCheck_ratio");
  //     Draw_TH1_Histogram(H1D_jetPt_ratio_mcpFoldedUnfoldedMcp, textContext, pdfName_ratio_McpFoldedUnfoldedCheck, texPtX, texRatioMcpFoldedUnfoldedMcp, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacementAuto, contextPlacementAuto, "ratioLine");
  //   }
  // }
}

// rename refoldedUnfolded as closure test?
// and try and spend 15 min to clean hist names for the spectrum analysis



// WARNING FOR EFFICIENCIES I SHOULD REREAD THIS BELOW!!
// hMcEfficiency_vsPt->Divide(hMcSignalCount_vsPt,TrueV0PtSpectrum_AnalysisBins, 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency


#endif