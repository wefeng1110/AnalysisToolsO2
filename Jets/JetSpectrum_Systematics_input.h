#ifndef JETSPECTRUM_SYSTEMATICS_INPUT_H
#define JETSPECTRUM_SYSTEMATICS_INPUT_H

#include "TFile.h"
#include "TString.h"
#include "JetSpectrum_inputs.h"

// 0-10% Pb-Pb systematics inputs. Reuse nDatasets and nominal inputs above.

//////////////////// Background fluctuations ////////////////////
TFile* file_O2Analysis_Bkg_Systematics_list[nDatasets] = {
  new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23zzh_bkg_DelatRs_HY683074/AnalysisResults.root")
};
const TString analysisWorkflowBkg_Systematics = "jet-background-analysis-task";
const TString analysisWorkflowBkg_DeltaR08 = "jet-background-analysis-task_DeltaR08";
const TString analysisWorkflowBkg_DeltaR15 = "jet-background-analysis-task_DeltaR15";

// Older background input, retained as a reference.
// TFile* file_O2Analysis_Bkg_Systematics_old_list[nDatasets] = {
//   new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23zzh_bkgAnalysis_sys_HY659248/AnalysisResults_25kHz.root")
// };

//////////////////// Tracking efficiency ////////////////////
TFile* file_O2Analysis_MCfile_TrackingEfficiency_Systematics[nDatasets] = {
  new TFile("/Users/fengwenhui/run_hyperloop/pbpb_jetjetMC/input/JJMC_cent0010_new_TrackingEff_sys_HY688109/AnalysisResults.root")
};
const TString analysisWorkflowMC_TrackingEfficiency_Systematics = "jet-spectra-charged_lead5";

// Older efficiency prescription, retained as a reference (not an independent source).
// TFile* file_O2Analysis_MCfile_TrackingEfficiency_Systematics_old[nDatasets] = {
//   new TFile("/Users/fengwenhui/run_hyperloop/pbpb_jetjetMC/input/JJMC_cent0010_trackingEfficiencySys_2timesRun2_HY684497/AnalysisResults.root")
// };

//////////////////// Event selection: |zVertex| < 8 cm ////////////////////
TFile* file_O2Analysis_ZVertex8_Systematics_list[nDatasets] = {
  new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23PbPbGolden_zVertex8_cent0010_1030_3050_HY669452/AnalysisResults_25kHz.root")
};
const TString analysisWorkflowData_ZVertex8_Systematics = "jet-spectra-charged_cent0010";
// No 0-10% zVertex-varied MC file was identified in the May inputs.

//////////////////// Event selection: centrality FT0C Variant1 ////////////////////
TFile* file_O2Analysis_Centrality_Systematics_list[nDatasets] = {
  new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23PbPb_gloden_centFT0C_systematics_HY687336/AnalysisResults_centFT0CVariant1_25kHz_train687336.root")
};
const TString analysisWorkflowData_Centrality_Systematics = "jet-spectra-charged_cent0010";

 
#endif
