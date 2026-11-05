#if defined(__CLING__)
R__LOAD_LIBRARY(libRooUnfold)
#endif

#include "TStyle.h"
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
#include <RooUnfold.h>
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

#include<array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <sstream>
#include <string.h>
using namespace std;

// Misc utilities
void SetStyle_Systematics(Bool_t graypalette=kFALSE);
void LoadLibs_Systematics();



void Get_systematics_UnfoldMethod(TH1D* &uncertainty, TH1D* &preBarlow, int iDataset, int iRadius, const char* const* methods, const int* parameters, int nMethods, std::string options, TH1D** reference = nullptr);
void Draw_Systematics_UnfoldMethod(int iDataset, int iRadius, const char* const* methods, const int* parameters, int nMethods, std::string options);
void Get_systematics_SvdParameter(TH1D* &uncertainty, TH1D* &preBarlow, int iDataset, int iRadius, int nominalK, std::string options, TH1D** reference = nullptr);
void Draw_Systematics_SvdParameter(int iDataset, int iRadius, int nominalK, std::string options);

/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void JetSpectrum_systematics() {
  cout << "Apply Barlow test: " << (applyBarlowTest ? "true" : "false") << endl;
  // Load necessary libraries
  LoadLibs_Systematics();
  // Set the default style
  SetStyle_Systematics();

  // TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  // TString* Extra = new TString("");

  // gathers the analysis options in a single char[]

  int iDataset = 0;
  int iRadius = 0;


  char optionsAnalysis_withoutUnfoldingMethod[100] = "";
  snprintf(optionsAnalysis_withoutUnfoldingMethod, sizeof(optionsAnalysis_withoutUnfoldingMethod), "%s", unfoldingPrior);


  const int nUnfoldingMethods = 2;
  const char* unfoldingMethodList[nUnfoldingMethods] = {"Svd", "Bayes"}; // default is the first one in this list
  int unfoldParameterInputList[2] = {11, 7};

  Draw_Systematics_UnfoldMethod(iDataset, iRadius, unfoldingMethodList, unfoldParameterInputList, nUnfoldingMethods, optionsAnalysis_withoutUnfoldingMethod);
  Draw_Systematics_SvdParameter(iDataset, iRadius, unfoldParameterInputList[0], optionsAnalysis_withoutUnfoldingMethod);
}

/////////////////////////////////////////////////////
/////////////////// Misc utilities //////////////////
/////////////////////////////////////////////////////

void LoadLibs_Systematics() {
  // gSystem->Load("libCore.so");  
  // gSystem->Load("libGeom.so");
  // gSystem->Load("libPhysics.so");
  // gSystem->Load("libVMC");
  // gSystem->Load("libTree");
  // gSystem->Load("libMinuit");
  // gSystem->Load("libSTEERBase");
  // gSystem->Load("libESD");
  // gSystem->Load("libAOD");
  // gSystem->Load("libANALYSIS");
  // gSystem->Load("libANALYSISalice");
  // gSystem->Load("libCORRFW");
  // gSystem->Load("libPWGTools");
}

void SetStyle_Systematics(Bool_t graypalette) {
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


// First spectrum is nominal. Use the largest absolute variation in each bin.
// When Barlow is enabled, test each variation before taking the accepted envelope.
void Get_Systematics_Envelope(TH1D* &uncertainty, TH1D* &preBarlow,
                              const std::vector<TH1D*>& spectra, TString name) {
  uncertainty = (TH1D*)spectra[0]->Clone("hSystematicUncertainty_" + name);
  uncertainty->Reset("ICES");
  preBarlow = nullptr;
  if (applyBarlowTest) {
    preBarlow = (TH1D*)uncertainty->Clone("hSystematicUncertainty_PreBarlow_" + name);
  }

  for (int bin = 1; bin <= spectra[0]->GetNbinsX(); ++bin) {
    double maximum = 0., acceptedMaximum = 0.;
    double maximumSigma = 0., acceptedSigma = 0.;
    for (size_t variation = 1; variation < spectra.size(); ++variation) {
      double difference = std::abs(spectra[variation]->GetBinContent(bin)
                                   - spectra[0]->GetBinContent(bin));
      if (!std::isfinite(difference)) throw std::runtime_error("Non-finite systematic variation");
      if (!applyBarlowTest) {
        maximum = std::max(maximum, difference);
        continue; // No Barlow error calculation or threshold test.
      }
      double error = spectra[variation]->GetBinError(bin);
      double referenceError = spectra[0]->GetBinError(bin);
      double sigma = std::sqrt(std::abs(error*error - referenceError*referenceError));
      if (!std::isfinite(sigma)) throw std::runtime_error("Non-finite Barlow uncertainty");
      if (difference > maximum) {
        maximum = difference;
        maximumSigma = sigma;
      }
      if (difference > N_SigmaBarlow*sigma && difference > acceptedMaximum) {
        acceptedMaximum = difference;
        acceptedSigma = sigma;
      }
    }
    uncertainty->SetBinContent(bin, applyBarlowTest ? acceptedMaximum : maximum);
    uncertainty->SetBinError(bin, applyBarlowTest ? acceptedSigma : 0.);
    if (applyBarlowTest) {
      preBarlow->SetBinContent(bin, maximum);
      preBarlow->SetBinError(bin, maximumSigma);
    }
  }
}

TH1D* Get_Systematics_MeasuredInput(int iDataset, int iRadius, std::string options) {
  TH1D* measured = nullptr;
  if (!normGenAndMeasByNEvtsForUnfoldingInput) {
    if (useFineBinningTest)
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measured, iDataset, iRadius, options);
    else
      Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measured, iDataset, iRadius, options);
  } else {
    if (useFineBinningTest)
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measured, iDataset, iRadius, options);
    else
      Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measured, iDataset, iRadius, options);
  }
  return measured;
}

void Draw_Systematics_Spectra(std::vector<TH1D*>& spectra, const std::vector<TString>& legends,
                             TString name, int iRadius) {
  TString* pdfName = new TString(name);
  TString* yTitle = normaliseUnfoldingResultsAtEnd
    ? texJet_d2Ndptdeta_EventNorm
    : new TString("d^{2}#it{N}_{jet}/d#it{p}_{T}d#it{#eta} (GeV/#it{c})^{-1}");
  std::array<std::array<float, 2>, 2> legendPlacement = {{{0.53, 0.59}, {0.94, 0.79}}};
  Draw_TH1_Histograms(spectra.data(), legends.data(), spectra.size(),
                     contextJetRadius(arrayRadius[iRadius]), pdfName, texPtJetGen, yTitle,
                     texCollisionDataInfo, drawnWindowUnfoldedMeasurement,
                     legendPlacement, contextPlacementAuto, "logy");
}

// Convert the absolute envelope to a fraction of the nominal spectrum.
// No independent ratio-error propagation: nominal and variations share the data.
void Draw_Systematics_RelativeUncertainty(TH1D* uncertainty, TH1D* preBarlow,
                                         TH1D* reference, TString name, int iRadius) {
  TH1D* histograms[2] = {uncertainty, preBarlow};
  int nHistograms = applyBarlowTest ? 2 : 1;
  for (int h = 0; h < nHistograms; ++h) {
    for (int bin = 1; bin <= reference->GetNbinsX(); ++bin) {
      double denominator = std::abs(reference->GetBinContent(bin));
      if (!std::isfinite(denominator)) throw std::runtime_error("Non-finite nominal spectrum");
      if (denominator == 0.) {
        if (histograms[h]->GetBinContent(bin) != 0.)
          throw std::runtime_error("Cannot form relative uncertainty with zero nominal yield");
        histograms[h]->SetBinContent(bin, 0.);
        histograms[h]->SetBinError(bin, 0.);
      } else {
        histograms[h]->SetBinContent(bin, histograms[h]->GetBinContent(bin)/denominator);
        histograms[h]->SetBinError(bin, histograms[h]->GetBinError(bin)/denominator);
      }
    }
    TString suffix = !applyBarlowTest ? "_NoBarlow" : (h == 0 ? "_Barlow" : "_PreBarlow");
    TString* pdfName = new TString(name + suffix);
    TString* yTitle = new TString("Relative systematic uncertainty");
    Draw_TH1_Histogram(histograms[h], contextJetRadius(arrayRadius[iRadius]), pdfName,
                       texPtJetGen, yTitle, texCollisionDataInfo, drawnWindowUnfoldedMeasurement,
                       legendPlacementAuto, contextPlacementAuto, "");
    cout << name << suffix << ": relative maximum = " << histograms[h]->GetMaximum() << endl;
  }
}

void Get_systematics_UnfoldMethod(TH1D* &uncertainty, TH1D* &preBarlow, int iDataset, int iRadius,
                                  const char* const* methods, const int* parameters, int nMethods,
                                  std::string options, TH1D** reference) {
  if (nMethods < 2) throw std::runtime_error("At least two unfolding methods are required");
  TH1D* measured = Get_Systematics_MeasuredInput(iDataset, iRadius, options);
  std::vector<TH1D*> spectra(nMethods);
  std::vector<TString> legends(nMethods);
  TString name = "UnfoldMethod_" + Datasets[iDataset] + Form("_R=%.1f", arrayRadius[iRadius]);
  for (int method = 0; method < nMethods; ++method) {
    Get_Pt_spectrum_unfolded(spectra[method], measured, iDataset, iRadius, parameters[method],
                             options + "," + methods[method]);
    legends[method] = TString(methods[method]) == "Svd"
      ? Form("SVD, k = %d", parameters[method])
      : Form("%s, iter. = %d", methods[method], parameters[method]);
    name += Form("_%s%d", methods[method], parameters[method]);
  }
  Get_Systematics_Envelope(uncertainty, preBarlow, spectra, name);
  Draw_Systematics_Spectra(spectra, legends, "Systematics_Spectra_" + name, iRadius);
  if (reference) *reference = spectra[0];
}

void Draw_Systematics_UnfoldMethod(int iDataset, int iRadius, const char* const* methods,
                                   const int* parameters, int nMethods, std::string options) {
  TH1D *uncertainty, *preBarlow, *reference;
  Get_systematics_UnfoldMethod(uncertainty, preBarlow, iDataset, iRadius, methods, parameters,
                               nMethods, options, &reference);
  TString name = "Systematics_UnfoldMethod_" + Datasets[iDataset]
                 + Form("_R=%.1f", arrayRadius[iRadius]);
  for (int method = 0; method < nMethods; ++method)
    name += Form("_%s%d", methods[method], parameters[method]);
  Draw_Systematics_RelativeUncertainty(uncertainty, preBarlow, reference, name, iRadius);
}

void Get_systematics_SvdParameter(TH1D* &uncertainty, TH1D* &preBarlow, int iDataset,
                                 int iRadius, int nominalK, std::string options, TH1D** reference) {
  if (automaticBestSvdParameter)
    throw std::runtime_error("Disable automaticBestSvdParameter for an explicit k variation");
  int nBins = useFineBinningTest ? nBinPtJetsFine[iRadius]
                               : std::min(nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius]);
  if (nominalK < 2 || nominalK + 1 > nBins)
    throw std::runtime_error("SVD k-1 and k+1 must lie within the response bin count");
  TH1D* measured = Get_Systematics_MeasuredInput(iDataset, iRadius, options);
  const int parameters[3] = {nominalK, nominalK - 1, nominalK + 1};
  std::vector<TH1D*> spectra(3);
  std::vector<TString> legends(3);
  TString name = "SvdParameter_" + Datasets[iDataset]
                 + Form("_R=%.1f_k%d", arrayRadius[iRadius], nominalK);
  for (int variation = 0; variation < 3; ++variation) {
    Get_Pt_spectrum_unfolded(spectra[variation], measured, iDataset, iRadius, parameters[variation],
                             options + ",Svd");
    legends[variation] = Form("SVD, k = %d%s", parameters[variation], variation == 0 ? " (nominal)" : "");
  }
  Get_Systematics_Envelope(uncertainty, preBarlow, spectra, name);
  Draw_Systematics_Spectra(spectra, legends, "Systematics_Spectra_" + name, iRadius);
  if (reference) *reference = spectra[0];
}

void Draw_Systematics_SvdParameter(int iDataset, int iRadius, int nominalK, std::string options) {
  TH1D *uncertainty, *preBarlow, *reference;
  Get_systematics_SvdParameter(uncertainty, preBarlow, iDataset, iRadius, nominalK, options, &reference);
  TString name = "Systematics_SvdParameter_" + Datasets[iDataset]
                 + Form("_R=%.1f_k%d", arrayRadius[iRadius], nominalK);
  Draw_Systematics_RelativeUncertainty(uncertainty, preBarlow, reference, name, iRadius);
}
