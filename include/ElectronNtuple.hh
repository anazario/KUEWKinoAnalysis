#ifndef ElectronNtuple_h
#define ElectronNtuple_h

#include "NtupleBase.hh"
#include "SUSYNANOBase.hh"
//#include "RestFrames/RestFrames.hh"

template class std::vector<std::vector<int> >;

//using namespace RestFrames;

class ElectronNtuple : public NtupleBase<SUSYNANOBase> {

public:
  ElectronNtuple(TTree* tree = 0);
  virtual ~ElectronNtuple();

private:
  TTree* InitOutputTree(const string& sample);
  void FillOutputTree(TTree* tree);

  void ClearVariables();

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

  int m_genNele;

  vector<double> m_genPT_ele;
  vector<double> m_genEta_ele;
  vector<double> m_genPhi_ele;
  vector<double> m_genM_ele;
  vector<int>    m_genCharge_ele;
  vector<int>    m_genPDGID_ele;
  vector<int>    m_genMomPDGID_ele;
  vector<int>    m_genIndex_ele;

  //nanoSUSY variables
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
  vector<bool>    m_Electron_mvaFall17V1Iso_WP80;
  vector<bool>    m_Electron_mvaFall17V1Iso_WP90;
  vector<bool>    m_Electron_mvaFall17V1Iso_WPL; 
  vector<bool>    m_Electron_mvaFall17V1noIso_WP80;
  vector<bool>    m_Electron_mvaFall17V1noIso_WP90;
  vector<bool>    m_Electron_mvaFall17V1noIso_WPL; 
  vector<bool>    m_Electron_mvaFall17V2Iso_WP80;  
  vector<bool>    m_Electron_mvaFall17V2Iso_WP90;  
  vector<bool>    m_Electron_mvaFall17V2Iso_WPL;   
  vector<bool>    m_Electron_mvaFall17V2noIso_WP80;
  vector<bool>    m_Electron_mvaFall17V2noIso_WP90;
  vector<bool>    m_Electron_mvaFall17V2noIso_WPL; 
  vector<bool>    m_Electron_mvaSpring16GP_WP80;   
  vector<bool>    m_Electron_mvaSpring16GP_WP90;   
  vector<bool>    m_Electron_mvaSpring16HZZ_WPL;   

  vector<int>     m_Electron_genPartIdx;   
  vector<char>    m_Electron_genPartFlav;  

  vector<UChar_t> m_Electron_cleanmask;   

};

#endif
