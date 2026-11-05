#ifndef JETSPECTRU_RESPONSEMATRIXFUNCTIONS_H
#define JETSPECTRU_RESPONSEMATRIXFUNCTIONS_H

void Get_PtResponseMatrix_Fluctuations(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int randomConeTypeInput, int iDataset, int iRadius, std::string optionsFluctResp);
void Get_PtResponseMatrix_Fluctuations(TH2D* &H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, std::string optionsFluctResp);
void Get_PtResponseMatrix_detectorResponse(TH2D* &H2D_jetPtResponseMatrix_detectorResponse, int iDataset, int iRadius, std::string options, bool controlMC);
void Get_PtResponseMatrix_DetectorAndFluctuationsCombined_postFinalise(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, std::string options, bool controlMC);
void Get_PtResponseMatrix_DetectorAndFluctuationsCombined_preFinalise(TH2D* &H2D_jetPtResponseMatrix_detectorAndFluctuationsCombined_fineBinning, TH2D* H2D_jetPtResponseMatrix_detectorResponse, TH2D* H2D_jetPtResponseMatrix_fluctuations, int iDataset, int iRadius, std::string options);
void ReweightResponseMatrixWithPrior(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options, bool controlMC);
void ReweightResponseMatrixWithPrior_fineBinningOnly(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options, bool controlMC);
void NormYSlicesAndscaleRespByXYWidth(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options);
void MergeResponseMatrixBins(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options);
void FinaliseResponseMatrix_priorAndNormYslicesAndMergeBins(TH2D* &H2D_jetPtResponseMatrix, int iDataset, int iRadius, std::string options, bool controlMC);

#endif