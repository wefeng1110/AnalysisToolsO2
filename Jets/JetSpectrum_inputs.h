
// This is a template. To use the JetSpectrum_DrawingMacro.C, rename this file to JetSpectrum_inputs.h and edit it how you want.


#ifndef JETSPECTRUM_INPUTS_H
#define JETSPECTRUM_INPUTS_H

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////       file access choice       ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
TFile* file_O2Analysis_run2ComparisonFileHannaBossiLaura = new TFile("/Users/fengwenhui/workspace/Unfolding_AreaBased_HannahMethod_R020_Nominal_ExtendedPtRange.root");
TFile* file_O2Analysis_run2ComparisonFileMLPaper = new TFile("/Users/fengwenhui/workspace/Ch-jetSuppression_PbPb502TeV.root");


//////// -------- LHC23zzh pass 4 with - pp sim anchored to PbPb 10% - lead05 ///////
TString* texEnergyPbPb = new TString("#sqrt{#it{s}_{NN}} = 5.36 TeV"); 
TString* texEnergy = new TString("pp, #sqrt{#it{s}} = 5.36 TeV"); 
TString* texCollisionDataType = new TString("0#font[122]{-}10% Pb#font[122]{-}Pb"); 
TString* texCollisionDataInfo = new TString((TString)*texCollisionDataType+", "+(TString)*texEnergyPbPb); 
TString* texCollisionMCType = new TString("PYTHIA + GEANT4"); 
TString* texCollisionMCInfo = new TString((TString)*texCollisionMCType+", "+(TString)*texEnergy); 
const TString* texDatasetsComparisonType = new TString("0#font[122]{-}10% cent.");
const TString* texDatasetsComparisonCommonDenominator = new TString("ALICE performance");

const int nDatasets = 1;
const TString Datasets[nDatasets] = {"LHC23PbPb_golden"
                                    };
const TString DatasetsNames[nDatasets] = {"Run3PbPb"};
TFile* file_O2Analysis_list[nDatasets] = {new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23pbpbGolden_25kHz_lead3_5_HY669451/AnalysisResults_jetanalysis_25kHzMax_goldenRuns_0010cent_nominal.root")
                                      };
// Nominal background fluctuations for the same 0-10% data selection.
TFile* file_O2Analysis_Bkg_list[nDatasets] = {new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23pbpbGolden_25kHz_lead3_5_HY669451/AnalysisResults_bkganalysis_25kHzMax_goldenRuns_0010cent_1030cent_5070cent.root")};
TFile* file_O2Analysis_MCfile_GeneralResponse[nDatasets] = {new TFile("/Users/fengwenhui/run_hyperloop/pbpb_jetjetMC/input/JJMC_good_cent0010_TrackEffReduction_HY685161/AnalysisResults.root")
                                                    };
TFile* file_O2Analysis_MCfile_UnfoldingControl_input[nDatasets] = {new TFile("/Users/fengwenhui/run_hyperloop/pbpb_jetjetMC/input/JJMC_good_cent0010_TrackEffReduction_HY685161/AnalysisResults.root")
                                                                        }; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)

TFile* file_O2Analysis_MCfile_UnfoldingControl_response[nDatasets] = {new TFile("/Users/fengwenhui/run_hyperloop/pbpb_jetjetMC/input/JJMC_good_cent0010_TrackEffReduction_HY685161/AnalysisResults.root")
                                                                        }; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)


const TString trainIdData = "";
const TString analysisWorkflowData = "jet-spectra-charged_cent0010"+trainIdData;
const bool collHistIsObsolete[nDatasets] = {false};

const TString trainIdBkg = "";
const TString analysisWorkflowBkg = "jet-background-analysis-task"+trainIdBkg;
const TString trainIdUnfoldingControl = "";
const TString analysisWorkflow_unfoldingControl = "jet-spectra-charged_lead5"+trainIdUnfoldingControl;

const TString trainIdMC = "";
const TString analysisWorkflowMC = "jet-spectra-charged_lead5"+trainIdMC;
// May 2026 MC stores matchedgeo_mcpetaconstraint (not mcpmcdetaconstraint).
const bool etaCutOnMatchedJetsIsObsoleteVersion = true;
const bool mcCollHistIsObsolete = false;


#endif
