#ifndef LeptonNtuple_h
#define LeptonNtuple_h

#include "NtupleBase.hh"
#include "SUSYNANOBase.hh"
//#include "RestFrames/RestFrames.hh"

template class std::vector<std::vector<int> >;

//using namespace RestFrames;

class LeptonNtuple : public NtupleBase<SUSYNANOBase> {

public:
  LeptonNtuple(TTree* tree = 0);
  virtual ~LeptonNtuple();

private:
  TTree* InitOutputTree(const string& sample);
  void FillOutputTree(TTree* tree);

  void ClearVariables();
  void InitNanoSUSYBranches();

  // common variables for output tree
  double m_weight;
  
  int m_Nele;
  
  vector<double> m_PT_ele;
  vector<double> m_Eta_ele;
  vector<double> m_Phi_ele;
  vector<double> m_M_ele;
  vector<int>    m_Charge_ele;
  vector<int>    m_PDGID_ele;
  vector<double> m_RelIso_ele;
  vector<double> m_MiniIso_ele;
  vector<int>    m_ID_ele;
  vector<int>    m_Index_ele;
  vector<bool>   m_baseline_loose_ele;
  vector<bool>   m_baseline_tight_ele;
  vector<bool>   m_baseline_ele;

  int m_Nmu;

  vector<double> m_PT_mu;
  vector<double> m_Eta_mu;
  vector<double> m_Phi_mu;
  vector<double> m_M_mu;
  vector<int>    m_Charge_mu;
  vector<int>    m_PDGID_mu;
  vector<double> m_RelIso_mu;
  vector<double> m_MiniIso_mu;
  vector<int>    m_ID_mu;
  vector<int>    m_Index_mu;
  vector<bool>   m_baseline_loose_mu;
  vector<bool>   m_baseline_tight_mu;

  int m_genNele;

  vector<double> m_genPT_ele;
  vector<double> m_genEta_ele;
  vector<double> m_genPhi_ele;
  vector<double> m_genM_ele;
  vector<int>    m_genCharge_ele;
  vector<int>    m_genPDGID_ele;
  vector<int>    m_genMomPDGID_ele;
  vector<int>    m_genIndex_ele;

  int m_genNmu;

  vector<double> m_genPT_mu;
  vector<double> m_genEta_mu;
  vector<double> m_genPhi_mu;
  vector<double> m_genM_mu;
  vector<int>    m_genCharge_mu;
  vector<int>    m_genPDGID_mu;
  vector<int>    m_genMomPDGID_mu;
  vector<int>    m_genIndex_mu;

  //nanoSUSY electron variables
  vector<float>   m_Electron_deltaEtaSC;
  vector<float>   m_Electron_dr03EcalRecHitSumEt;
  vector<float>   m_Electron_dr03HcalDepth1TowerSumEt;
  vector<float>   m_Electron_dr03TkSumPt;
  vector<float>   m_Electron_dr03TkSumPtHEEP;
  vector<float>   m_Electron_dxy;
  vector<float>   m_Electron_dxyErr;
  vector<float>   m_Electron_dz;  
  vector<float>   m_Electron_dzErr;
  vector<float>   m_Electron_eCorr;
  vector<float>   m_Electron_eInvMinusPInv;
  vector<float>   m_Electron_energyErr;   
  vector<float>   m_Electron_eta;   
  vector<float>   m_Electron_hoe;   
  vector<float>   m_Electron_ip3d;  
  vector<float>   m_Electron_jetRelIso;
  vector<float>   m_Electron_mass;   
  vector<float>   m_Electron_miniPFRelIso_all;
  vector<float>   m_Electron_miniPFRelIso_chg;
  vector<float>   m_Electron_mvaFall17V1Iso;  
  vector<float>   m_Electron_mvaFall17V1noIso;
  vector<float>   m_Electron_mvaFall17V2Iso;  
  vector<float>   m_Electron_mvaFall17V2noIso;
  vector<float>   m_Electron_mvaSpring16GP;   
  vector<float>   m_Electron_mvaSpring16HZZ;  
  vector<float>   m_Electron_pfRelIso03_all;  
  vector<float>   m_Electron_pfRelIso03_chg;  
  vector<float>   m_Electron_phi;  
  vector<float>   m_Electron_pt;   
  vector<float>   m_Electron_r9;   
  vector<float>   m_Electron_sieie;
  vector<float>   m_Electron_sip3d;
  vector<float>   m_Electron_mvaTTH;
  vector<int>     m_Electron_charge;
  vector<int>     m_Electron_cutBased;
  vector<int>     m_Electron_cutBased_Fall17_V1;
  vector<int>     m_Electron_cutBased_HLTPreSel;
  vector<int>     m_Electron_cutBased_Spring15; 
  vector<int>     m_Electron_cutBased_Sum16;
  vector<int>     m_Electron_jetIdx;
  vector<int>     m_Electron_pdgId; 
  vector<int>     m_Electron_photonIdx; 
  vector<int>     m_Electron_tightCharge;
  vector<int>     m_Electron_vidNestedWPBitmap;
  vector<int>     m_Electron_vidNestedWPBitmapSpring15; 
  vector<int>     m_Electron_vidNestedWPBitmapSum16;   
  vector<bool>    m_Electron_convVeto; 
  vector<bool>    m_Electron_cutBased_HEEP; 
  vector<bool>    m_Electron_isPFcand;
  vector<char>    m_Electron_lostHits;
  vector<bool>    m_Electron_mvaFall17V1noIso_WP80;
  vector<bool>    m_Electron_mvaFall17V1noIso_WP90;
  vector<bool>    m_Electron_mvaFall17V1noIso_WPL; 
  vector<bool>    m_Electron_mvaFall17V2noIso_WP80;
  vector<bool>    m_Electron_mvaFall17V2noIso_WP90;
  vector<bool>    m_Electron_mvaFall17V2noIso_WPL; 
  vector<bool>    m_Electron_mvaSpring16GP_WP80;   
  vector<bool>    m_Electron_mvaSpring16GP_WP90;   
  vector<bool>    m_Electron_mvaSpring16HZZ_WPL;   

  vector<int>     m_Electron_genPartIdx;   
  vector<char>    m_Electron_genPartFlav;  

  //nanoSUSY muon variables
  vector<float>   m_Muon_dxy;
  vector<float>   m_Muon_dxyErr;
  vector<float>   m_Muon_dz;
  vector<float>   m_Muon_dzErr;
  vector<float>   m_Muon_eta;
  vector<float>   m_Muon_ip3d;
  vector<float>   m_Muon_jetRelIso;
  vector<float>   m_Muon_mass;
  vector<float>   m_Muon_miniPFRelIso_all;
  vector<float>   m_Muon_miniPFRelIso_chg;
  vector<float>   m_Muon_pfRelIso03_all;
  vector<float>   m_Muon_pfRelIso03_chg;
  vector<float>   m_Muon_pfRelIso04_all;
  vector<float>   m_Muon_phi;
  vector<float>   m_Muon_pt;
  vector<float>   m_Muon_ptErr;
  vector<float>   m_Muon_segmentComp;
  vector<float>   m_Muon_sip3d;
  vector<float>   m_Muon_mvaTTH;
  vector<int>     m_Muon_charge;
  vector<int>     m_Muon_jetIdx;
  vector<int>     m_Muon_nStations;
  vector<int>     m_Muon_nTrackerLayers;
  vector<int>     m_Muon_pdgId;
  vector<int>     m_Muon_tightCharge;
  vector<char>    m_Muon_highPtId;
  vector<bool>    m_Muon_inTimeMuon;
  vector<bool>    m_Muon_isGlobal;
  vector<bool>    m_Muon_isPFcand;
  vector<bool>    m_Muon_isTracker;
  vector<bool>    m_Muon_mediumId;
  vector<bool>    m_Muon_mediumPromptId;
  vector<bool>    m_Muon_miniIsoId;
  vector<bool>    m_Muon_multiIsoId;
  vector<bool>    m_Muon_mvaId;
  vector<char>    m_Muon_pfIsoId;
  vector<bool>    m_Muon_softId;
  vector<bool>    m_Muon_softMvaId;
  vector<bool>    m_Muon_tightId;
  vector<char>    m_Muon_tkIsoId;
  vector<bool>    m_Muon_triggerIdLoose;

};

#endif
