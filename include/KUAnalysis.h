//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Dec  8 15:33:37 2019 by ROOT version 6.14/09
// from TTree KUAnalysis/KUAnalysis
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef KUAnalysis_h
#define KUAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;

class KUAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        weight;
   Double_t        MET_pt;
   Int_t           nGenPart;
   vector<float>   *GenPart_eta;
   vector<float>   *GenPart_mass;
   vector<float>   *GenPart_phi;
   vector<float>   *GenPart_pt;
   vector<int>     *GenPart_genPartIdxMother;
   vector<int>     *GenPart_pdgId;
   vector<int>     *GenPart_status;
   Int_t           Nele;
   vector<double>  *PT_ele;
   vector<double>  *Eta_ele;
   vector<double>  *Phi_ele;
   vector<double>  *M_ele;
   vector<int>     *Charge_ele;
   vector<int>     *PDGID_ele;
   vector<double>  *RelIso_ele;
   vector<double>  *MiniIso_ele;
   vector<int>     *ID_ele;
   vector<int>     *Index_ele;
   vector<bool>    *FO_baseline_ele;
   vector<bool>    *baseline_loose_ele;
   vector<bool>    *baseline_tight_ele;
   Int_t           genNele;
   vector<double>  *genPT_ele;
   vector<double>  *genEta_ele;
   vector<double>  *genPhi_ele;
   vector<double>  *genM_ele;
   vector<int>     *genCharge_ele;
   vector<int>     *genPDGID_ele;
   vector<int>     *genMomPDGID_ele;
   vector<int>     *genIndex_ele;
   Int_t           Nmu;
   vector<double>  *PT_mu;
   vector<double>  *Eta_mu;
   vector<double>  *Phi_mu;
   vector<double>  *M_mu;
   vector<int>     *Charge_mu;
   vector<int>     *PDGID_mu;
   vector<double>  *RelIso_mu;
   vector<double>  *MiniIso_mu;
   vector<int>     *ID_mu;
   vector<int>     *Index_mu;
   vector<bool>    *baseline_loose_mu;
   vector<bool>    *baseline_tight_mu;
   Int_t           genNmu;
   vector<double>  *genPT_mu;
   vector<double>  *genEta_mu;
   vector<double>  *genPhi_mu;
   vector<double>  *genM_mu;
   vector<int>     *genCharge_mu;
   vector<int>     *genPDGID_mu;
   vector<int>     *genMomPDGID_mu;
   vector<int>     *genIndex_mu;
   vector<float>   *Electron_deltaEtaSC;
   vector<float>   *Electron_dr03EcalRecHitSumEt;
   vector<float>   *Electron_dr03HcalDepth1TowerSumEt;
   vector<float>   *Electron_dr03TkSumPt;
   vector<float>   *Electron_dr03TkSumPtHEEP;
   vector<float>   *Electron_dxy;
   vector<float>   *Electron_dxyErr;
   vector<float>   *Electron_dz;
   vector<float>   *Electron_dzErr;
   vector<float>   *Electron_eCorr;
   vector<float>   *Electron_eInvMinusPInv;
   vector<float>   *Electron_energyErr;
   vector<float>   *Electron_eta;
   vector<float>   *Electron_hoe;
   vector<float>   *Electron_ip3d;
   vector<float>   *Electron_jetRelIso;
   vector<float>   *Electron_mass;
   vector<float>   *Electron_miniPFRelIso_all;
   vector<float>   *Electron_miniPFRelIso_chg;
   vector<float>   *Electron_mvaFall17V1Iso;
   vector<float>   *Electron_mvaFall17V1noIso;
   vector<float>   *Electron_mvaFall17V2Iso;
   vector<float>   *Electron_mvaFall17V2noIso;
   vector<float>   *Electron_pfRelIso03_all;
   vector<float>   *Electron_pfRelIso03_chg;
   vector<float>   *Electron_phi;
   vector<float>   *Electron_pt;
   vector<float>   *Electron_r9;
   vector<float>   *Electron_sieie;
   vector<float>   *Electron_sip3d;
   vector<float>   *Electron_mvaTTH;
   vector<int>     *Electron_charge;
   vector<int>     *Electron_cutBased;
   vector<int>     *Electron_cutBased_Fall17_V1;
   vector<int>     *Electron_jetIdx;
   vector<int>     *Electron_pdgId;
   vector<int>     *Electron_photonIdx;
   vector<int>     *Electron_tightCharge;
   vector<int>     *Electron_vidNestedWPBitmap;
   vector<bool>    *Electron_convVeto;
   vector<bool>    *Electron_cutBased_HEEP;
   vector<bool>    *Electron_isPFcand;
   vector<char>    *Electron_lostHits;
   vector<bool>    *Electron_mvaFall17V1noIso_WP80;
   vector<bool>    *Electron_mvaFall17V1noIso_WP90;
   vector<bool>    *Electron_mvaFall17V1noIso_WPL;
   vector<bool>    *Electron_mvaFall17V2noIso_WP80;
   vector<bool>    *Electron_mvaFall17V2noIso_WP90;
   vector<bool>    *Electron_mvaFall17V2noIso_WPL;
   vector<int>     *Electron_genPartIdx;
   vector<char>    *Electron_genPartFlav;
   vector<float>   *Muon_dxy;
   vector<float>   *Muon_dxyErr;
   vector<float>   *Muon_dz;
   vector<float>   *Muon_dzErr;
   vector<float>   *Muon_eta;
   vector<float>   *Muon_ip3d;
   vector<float>   *Muon_jetRelIso;
   vector<float>   *Muon_mass;
   vector<float>   *Muon_miniPFRelIso_all;
   vector<float>   *Muon_miniPFRelIso_chg;
   vector<float>   *Muon_pfRelIso03_all;
   vector<float>   *Muon_pfRelIso03_chg;
   vector<float>   *Muon_pfRelIso04_all;
   vector<float>   *Muon_phi;
   vector<float>   *Muon_pt;
   vector<float>   *Muon_ptErr;
   vector<float>   *Muon_segmentComp;
   vector<float>   *Muon_sip3d;
   vector<float>   *Muon_mvaTTH;
   vector<int>     *Muon_charge;
   vector<int>     *Muon_jetIdx;
   vector<int>     *Muon_nStations;
   vector<int>     *Muon_nTrackerLayers;
   vector<int>     *Muon_pdgId;
   vector<int>     *Muon_tightCharge;
   vector<char>    *Muon_highPtId;
   vector<bool>    *Muon_inTimeMuon;
   vector<bool>    *Muon_isGlobal;
   vector<bool>    *Muon_isPFcand;
   vector<bool>    *Muon_isTracker;
   vector<bool>    *Muon_mediumId;
   vector<bool>    *Muon_mediumPromptId;
   vector<bool>    *Muon_miniIsoId;
   vector<bool>    *Muon_multiIsoId;
   vector<bool>    *Muon_mvaId;
   vector<char>    *Muon_pfIsoId;
   vector<bool>    *Muon_softId;
   vector<bool>    *Muon_softMvaId;
   vector<bool>    *Muon_tightId;
   vector<char>    *Muon_tkIsoId;
   vector<bool>    *Muon_triggerIdLoose;

   // List of branches
   TBranch        *b_weight;   //!
   TBranch        *b_MET_pt;   //!
   TBranch        *b_nGenPart;   //!                                                                                                                
   TBranch        *b_GenPart_eta;   //!                                                                                                             
   TBranch        *b_GenPart_mass;   //!                                                                                                            
   TBranch        *b_GenPart_phi;   //!                                                                                                             
   TBranch        *b_GenPart_pt;   //!                                                                                                              
   TBranch        *b_GenPart_genPartIdxMother;   //!                                                                                                
   TBranch        *b_GenPart_pdgId;   //!                                                                                                           
   TBranch        *b_GenPart_status;   //!
   TBranch        *b_Nele;   //!
   TBranch        *b_PT_ele;   //!
   TBranch        *b_Eta_ele;   //!
   TBranch        *b_Phi_ele;   //!
   TBranch        *b_M_ele;   //!
   TBranch        *b_Charge_ele;   //!
   TBranch        *b_PDGID_ele;   //!
   TBranch        *b_RelIso_ele;   //!
   TBranch        *b_MiniIso_ele;   //!
   TBranch        *b_ID_ele;   //!
   TBranch        *b_Index_ele;   //!
   TBranch        *b_FO_baseline_ele;   //!
   TBranch        *b_baseline_loose_ele;   //!
   TBranch        *b_baseline_tight_ele;   //!
   TBranch        *b_genNele;   //!
   TBranch        *b_genPT_ele;   //!
   TBranch        *b_genEta_ele;   //!
   TBranch        *b_genPhi_ele;   //!
   TBranch        *b_genM_ele;   //!
   TBranch        *b_genCharge_ele;   //!
   TBranch        *b_genPDGID_ele;   //!
   TBranch        *b_genMomPDGID_ele;   //!
   TBranch        *b_genIndex_ele;   //!
   TBranch        *b_Nmu;   //!
   TBranch        *b_PT_mu;   //!
   TBranch        *b_Eta_mu;   //!
   TBranch        *b_Phi_mu;   //!
   TBranch        *b_M_mu;   //!
   TBranch        *b_Charge_mu;   //!
   TBranch        *b_PDGID_mu;   //!
   TBranch        *b_RelIso_mu;   //!
   TBranch        *b_MiniIso_mu;   //!
   TBranch        *b_ID_mu;   //!
   TBranch        *b_Index_mu;   //!
   TBranch        *b_baseline_loose_mu;   //!
   TBranch        *b_baseline_tight_mu;   //!
   TBranch        *b_genNmu;   //!
   TBranch        *b_genPT_mu;   //!
   TBranch        *b_genEta_mu;   //!
   TBranch        *b_genPhi_mu;   //!
   TBranch        *b_genM_mu;   //!
   TBranch        *b_genCharge_mu;   //!
   TBranch        *b_genPDGID_mu;   //!
   TBranch        *b_genMomPDGID_mu;   //!
   TBranch        *b_genIndex_mu;   //!
   TBranch        *b_Electron_deltaEtaSC;   //!
   TBranch        *b_Electron_dr03EcalRecHitSumEt;   //!
   TBranch        *b_Electron_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_Electron_dr03TkSumPt;   //!
   TBranch        *b_Electron_dr03TkSumPtHEEP;   //!
   TBranch        *b_Electron_dxy;   //!
   TBranch        *b_Electron_dxyErr;   //!
   TBranch        *b_Electron_dz;   //!
   TBranch        *b_Electron_dzErr;   //!
   TBranch        *b_Electron_eCorr;   //!
   TBranch        *b_Electron_eInvMinusPInv;   //!
   TBranch        *b_Electron_energyErr;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_hoe;   //!
   TBranch        *b_Electron_ip3d;   //!
   TBranch        *b_Electron_jetRelIso;   //!
   TBranch        *b_Electron_mass;   //!
   TBranch        *b_Electron_miniPFRelIso_all;   //!
   TBranch        *b_Electron_miniPFRelIso_chg;   //!
   TBranch        *b_Electron_mvaFall17V1Iso;   //!
   TBranch        *b_Electron_mvaFall17V1noIso;   //!
   TBranch        *b_Electron_mvaFall17V2Iso;   //!
   TBranch        *b_Electron_mvaFall17V2noIso;   //!
   TBranch        *b_Electron_pfRelIso03_all;   //!
   TBranch        *b_Electron_pfRelIso03_chg;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_r9;   //!
   TBranch        *b_Electron_sieie;   //!
   TBranch        *b_Electron_sip3d;   //!
   TBranch        *b_Electron_mvaTTH;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_cutBased;   //!
   TBranch        *b_Electron_cutBased_Fall17_V1;   //!
   TBranch        *b_Electron_jetIdx;   //!
   TBranch        *b_Electron_pdgId;   //!
   TBranch        *b_Electron_photonIdx;   //!
   TBranch        *b_Electron_tightCharge;   //!
   TBranch        *b_Electron_vidNestedWPBitmap;   //!
   TBranch        *b_Electron_convVeto;   //!
   TBranch        *b_Electron_cutBased_HEEP;   //!
   TBranch        *b_Electron_isPFcand;   //!
   TBranch        *b_Electron_lostHits;   //!
   TBranch        *b_Electron_mvaFall17V1noIso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V1noIso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V1noIso_WPL;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WPL;   //!
   TBranch        *b_Electron_genPartIdx;   //!
   TBranch        *b_Electron_genPartFlav;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_dxyErr;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_dzErr;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_ip3d;   //!
   TBranch        *b_Muon_jetRelIso;   //!
   TBranch        *b_Muon_mass;   //!
   TBranch        *b_Muon_miniPFRelIso_all;   //!
   TBranch        *b_Muon_miniPFRelIso_chg;   //!
   TBranch        *b_Muon_pfRelIso03_all;   //!
   TBranch        *b_Muon_pfRelIso03_chg;   //!
   TBranch        *b_Muon_pfRelIso04_all;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_ptErr;   //!
   TBranch        *b_Muon_segmentComp;   //!
   TBranch        *b_Muon_sip3d;   //!
   TBranch        *b_Muon_mvaTTH;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_jetIdx;   //!
   TBranch        *b_Muon_nStations;   //!
   TBranch        *b_Muon_nTrackerLayers;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_tightCharge;   //!
   TBranch        *b_Muon_highPtId;   //!
   TBranch        *b_Muon_inTimeMuon;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_isPFcand;   //!
   TBranch        *b_Muon_isTracker;   //!
   TBranch        *b_Muon_mediumId;   //!
   TBranch        *b_Muon_mediumPromptId;   //!
   TBranch        *b_Muon_miniIsoId;   //!
   TBranch        *b_Muon_multiIsoId;   //!
   TBranch        *b_Muon_mvaId;   //!
   TBranch        *b_Muon_pfIsoId;   //!
   TBranch        *b_Muon_softId;   //!
   TBranch        *b_Muon_softMvaId;   //!
   TBranch        *b_Muon_tightId;   //!
   TBranch        *b_Muon_tkIsoId;   //!
   TBranch        *b_Muon_triggerIdLoose;   //!

   KUAnalysis(TTree *tree=0);
   virtual ~KUAnalysis();

   virtual void     Init(TTree *tree);
};

#endif


KUAnalysis::KUAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("test.root");
      }
      f->GetObject("KUAnalysis",tree);

   }
   Init(tree);
}

KUAnalysis::~KUAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

inline void KUAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   GenPart_eta = 0;
   GenPart_mass = 0;
   GenPart_phi = 0;
   GenPart_pt = 0;
   GenPart_genPartIdxMother = 0;
   GenPart_pdgId = 0;
   GenPart_status = 0;
   PT_ele = 0;
   Eta_ele = 0;
   Phi_ele = 0;
   M_ele = 0;
   Charge_ele = 0;
   PDGID_ele = 0;
   RelIso_ele = 0;
   MiniIso_ele = 0;
   ID_ele = 0;
   Index_ele = 0;
   FO_baseline_ele = 0;
   baseline_loose_ele = 0;
   baseline_tight_ele = 0;
   genPT_ele = 0;
   genEta_ele = 0;
   genPhi_ele = 0;
   genM_ele = 0;
   genCharge_ele = 0;
   genPDGID_ele = 0;
   genMomPDGID_ele = 0;
   genIndex_ele = 0;
   PT_mu = 0;
   Eta_mu = 0;
   Phi_mu = 0;
   M_mu = 0;
   Charge_mu = 0;
   PDGID_mu = 0;
   RelIso_mu = 0;
   MiniIso_mu = 0;
   ID_mu = 0;
   Index_mu = 0;
   baseline_loose_mu = 0;
   baseline_tight_mu = 0;
   genPT_mu = 0;
   genEta_mu = 0;
   genPhi_mu = 0;
   genM_mu = 0;
   genCharge_mu = 0;
   genPDGID_mu = 0;
   genMomPDGID_mu = 0;
   genIndex_mu = 0;
   Electron_deltaEtaSC = 0;
   Electron_dr03EcalRecHitSumEt = 0;
   Electron_dr03HcalDepth1TowerSumEt = 0;
   Electron_dr03TkSumPt = 0;
   Electron_dr03TkSumPtHEEP = 0;
   Electron_dxy = 0;
   Electron_dxyErr = 0;
   Electron_dz = 0;
   Electron_dzErr = 0;
   Electron_eCorr = 0;
   Electron_eInvMinusPInv = 0;
   Electron_energyErr = 0;
   Electron_eta = 0;
   Electron_hoe = 0;
   Electron_ip3d = 0;
   Electron_jetRelIso = 0;
   Electron_mass = 0;
   Electron_miniPFRelIso_all = 0;
   Electron_miniPFRelIso_chg = 0;
   Electron_mvaFall17V1Iso = 0;
   Electron_mvaFall17V1noIso = 0;
   Electron_mvaFall17V2Iso = 0;
   Electron_mvaFall17V2noIso = 0;
   Electron_pfRelIso03_all = 0;
   Electron_pfRelIso03_chg = 0;
   Electron_phi = 0;
   Electron_pt = 0;
   Electron_r9 = 0;
   Electron_sieie = 0;
   Electron_sip3d = 0;
   Electron_mvaTTH = 0;
   Electron_charge = 0;
   Electron_cutBased = 0;
   Electron_cutBased_Fall17_V1 = 0;
   Electron_jetIdx = 0;
   Electron_pdgId = 0;
   Electron_photonIdx = 0;
   Electron_tightCharge = 0;
   Electron_vidNestedWPBitmap = 0;
   Electron_convVeto = 0;
   Electron_cutBased_HEEP = 0;
   Electron_isPFcand = 0;
   Electron_lostHits = 0;
   Electron_mvaFall17V1noIso_WP80 = 0;
   Electron_mvaFall17V1noIso_WP90 = 0;
   Electron_mvaFall17V1noIso_WPL = 0;
   Electron_mvaFall17V2noIso_WP80 = 0;
   Electron_mvaFall17V2noIso_WP90 = 0;
   Electron_mvaFall17V2noIso_WPL = 0;
   Electron_genPartIdx = 0;
   Electron_genPartFlav = 0;
   Muon_dxy = 0;
   Muon_dxyErr = 0;
   Muon_dz = 0;
   Muon_dzErr = 0;
   Muon_eta = 0;
   Muon_ip3d = 0;
   Muon_jetRelIso = 0;
   Muon_mass = 0;
   Muon_miniPFRelIso_all = 0;
   Muon_miniPFRelIso_chg = 0;
   Muon_pfRelIso03_all = 0;
   Muon_pfRelIso03_chg = 0;
   Muon_pfRelIso04_all = 0;
   Muon_phi = 0;
   Muon_pt = 0;
   Muon_ptErr = 0;
   Muon_segmentComp = 0;
   Muon_sip3d = 0;
   Muon_mvaTTH = 0;
   Muon_charge = 0;
   Muon_jetIdx = 0;
   Muon_nStations = 0;
   Muon_nTrackerLayers = 0;
   Muon_pdgId = 0;
   Muon_tightCharge = 0;
   Muon_highPtId = 0;
   Muon_inTimeMuon = 0;
   Muon_isGlobal = 0;
   Muon_isPFcand = 0;
   Muon_isTracker = 0;
   Muon_mediumId = 0;
   Muon_mediumPromptId = 0;
   Muon_miniIsoId = 0;
   Muon_multiIsoId = 0;
   Muon_mvaId = 0;
   Muon_pfIsoId = 0;
   Muon_softId = 0;
   Muon_softMvaId = 0;
   Muon_tightId = 0;
   Muon_tkIsoId = 0;
   Muon_triggerIdLoose = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("GenPart_eta", &GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_mass", &GenPart_mass, &b_GenPart_mass);
   fChain->SetBranchAddress("GenPart_phi", &GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_pt", &GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
   fChain->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId, &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_status", &GenPart_status, &b_GenPart_status);
   fChain->SetBranchAddress("Nele", &Nele, &b_Nele);
   fChain->SetBranchAddress("PT_ele", &PT_ele, &b_PT_ele);
   fChain->SetBranchAddress("Eta_ele", &Eta_ele, &b_Eta_ele);
   fChain->SetBranchAddress("Phi_ele", &Phi_ele, &b_Phi_ele);
   fChain->SetBranchAddress("M_ele", &M_ele, &b_M_ele);
   fChain->SetBranchAddress("Charge_ele", &Charge_ele, &b_Charge_ele);
   fChain->SetBranchAddress("PDGID_ele", &PDGID_ele, &b_PDGID_ele);
   fChain->SetBranchAddress("RelIso_ele", &RelIso_ele, &b_RelIso_ele);
   fChain->SetBranchAddress("MiniIso_ele", &MiniIso_ele, &b_MiniIso_ele);
   fChain->SetBranchAddress("ID_ele", &ID_ele, &b_ID_ele);
   fChain->SetBranchAddress("Index_ele", &Index_ele, &b_Index_ele);
   fChain->SetBranchAddress("FO_baseline_ele", &FO_baseline_ele, &b_FO_baseline_ele);
   fChain->SetBranchAddress("baseline_loose_ele", &baseline_loose_ele, &b_baseline_loose_ele);
   fChain->SetBranchAddress("baseline_tight_ele", &baseline_tight_ele, &b_baseline_tight_ele);
   fChain->SetBranchAddress("genNele", &genNele, &b_genNele);
   fChain->SetBranchAddress("genPT_ele", &genPT_ele, &b_genPT_ele);
   fChain->SetBranchAddress("genEta_ele", &genEta_ele, &b_genEta_ele);
   fChain->SetBranchAddress("genPhi_ele", &genPhi_ele, &b_genPhi_ele);
   fChain->SetBranchAddress("genM_ele", &genM_ele, &b_genM_ele);
   fChain->SetBranchAddress("genCharge_ele", &genCharge_ele, &b_genCharge_ele);
   fChain->SetBranchAddress("Electron_genPartIdx", &genPDGID_ele, &b_genPDGID_ele);
   fChain->SetBranchAddress("genMomPDGID_ele", &genMomPDGID_ele, &b_genMomPDGID_ele);
   fChain->SetBranchAddress("genIndex_ele", &genIndex_ele, &b_genIndex_ele);
   fChain->SetBranchAddress("Nmu", &Nmu, &b_Nmu);
   fChain->SetBranchAddress("PT_mu", &PT_mu, &b_PT_mu);
   fChain->SetBranchAddress("Eta_mu", &Eta_mu, &b_Eta_mu);
   fChain->SetBranchAddress("Phi_mu", &Phi_mu, &b_Phi_mu);
   fChain->SetBranchAddress("M_mu", &M_mu, &b_M_mu);
   fChain->SetBranchAddress("Charge_mu", &Charge_mu, &b_Charge_mu);
   fChain->SetBranchAddress("PDGID_mu", &PDGID_mu, &b_PDGID_mu);
   fChain->SetBranchAddress("RelIso_mu", &RelIso_mu, &b_RelIso_mu);
   fChain->SetBranchAddress("MiniIso_mu", &MiniIso_mu, &b_MiniIso_mu);
   fChain->SetBranchAddress("ID_mu", &ID_mu, &b_ID_mu);
   fChain->SetBranchAddress("Index_mu", &Index_mu, &b_Index_mu);
   fChain->SetBranchAddress("baseline_loose_mu", &baseline_loose_mu, &b_baseline_loose_mu);
   fChain->SetBranchAddress("baseline_tight_mu", &baseline_tight_mu, &b_baseline_tight_mu);
   fChain->SetBranchAddress("genNmu", &genNmu, &b_genNmu);
   fChain->SetBranchAddress("genPT_mu", &genPT_mu, &b_genPT_mu);
   fChain->SetBranchAddress("genEta_mu", &genEta_mu, &b_genEta_mu);
   fChain->SetBranchAddress("genPhi_mu", &genPhi_mu, &b_genPhi_mu);
   fChain->SetBranchAddress("genM_mu", &genM_mu, &b_genM_mu);
   fChain->SetBranchAddress("genCharge_mu", &genCharge_mu, &b_genCharge_mu);
   fChain->SetBranchAddress("genPDGID_mu", &genPDGID_mu, &b_genPDGID_mu);
   fChain->SetBranchAddress("genMomPDGID_mu", &genMomPDGID_mu, &b_genMomPDGID_mu);
   fChain->SetBranchAddress("genIndex_mu", &genIndex_mu, &b_genIndex_mu);
   fChain->SetBranchAddress("Electron_deltaEtaSC", &Electron_deltaEtaSC, &b_Electron_deltaEtaSC);
   fChain->SetBranchAddress("Electron_dr03EcalRecHitSumEt", &Electron_dr03EcalRecHitSumEt, &b_Electron_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", &Electron_dr03HcalDepth1TowerSumEt, &b_Electron_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Electron_dr03TkSumPt", &Electron_dr03TkSumPt, &b_Electron_dr03TkSumPt);
   fChain->SetBranchAddress("Electron_dr03TkSumPtHEEP", &Electron_dr03TkSumPtHEEP, &b_Electron_dr03TkSumPtHEEP);
   fChain->SetBranchAddress("Electron_dxy", &Electron_dxy, &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dxyErr", &Electron_dxyErr, &b_Electron_dxyErr);
   fChain->SetBranchAddress("Electron_dz", &Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron_dzErr", &Electron_dzErr, &b_Electron_dzErr);
   fChain->SetBranchAddress("Electron_eCorr", &Electron_eCorr, &b_Electron_eCorr);
   fChain->SetBranchAddress("Electron_eInvMinusPInv", &Electron_eInvMinusPInv, &b_Electron_eInvMinusPInv);
   fChain->SetBranchAddress("Electron_energyErr", &Electron_energyErr, &b_Electron_energyErr);
   fChain->SetBranchAddress("Electron_eta", &Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_hoe", &Electron_hoe, &b_Electron_hoe);
   fChain->SetBranchAddress("Electron_ip3d", &Electron_ip3d, &b_Electron_ip3d);
   fChain->SetBranchAddress("Electron_jetRelIso", &Electron_jetRelIso, &b_Electron_jetRelIso);
   fChain->SetBranchAddress("Electron_mass", &Electron_mass, &b_Electron_mass);
   fChain->SetBranchAddress("Electron_miniPFRelIso_all", &Electron_miniPFRelIso_all, &b_Electron_miniPFRelIso_all);
   fChain->SetBranchAddress("Electron_miniPFRelIso_chg", &Electron_miniPFRelIso_chg, &b_Electron_miniPFRelIso_chg);
   fChain->SetBranchAddress("Electron_mvaFall17V1Iso", &Electron_mvaFall17V1Iso, &b_Electron_mvaFall17V1Iso);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso", &Electron_mvaFall17V1noIso, &b_Electron_mvaFall17V1noIso);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso", &Electron_mvaFall17V2Iso, &b_Electron_mvaFall17V2Iso);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso", &Electron_mvaFall17V2noIso, &b_Electron_mvaFall17V2noIso);
   fChain->SetBranchAddress("Electron_pfRelIso03_all", &Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
   fChain->SetBranchAddress("Electron_pfRelIso03_chg", &Electron_pfRelIso03_chg, &b_Electron_pfRelIso03_chg);
   fChain->SetBranchAddress("Electron_phi", &Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", &Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_r9", &Electron_r9, &b_Electron_r9);
   fChain->SetBranchAddress("Electron_sieie", &Electron_sieie, &b_Electron_sieie);
   fChain->SetBranchAddress("Electron_sip3d", &Electron_sip3d, &b_Electron_sip3d);
   fChain->SetBranchAddress("Electron_mvaTTH", &Electron_mvaTTH, &b_Electron_mvaTTH);
   fChain->SetBranchAddress("Electron_charge", &Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_cutBased", &Electron_cutBased, &b_Electron_cutBased);
   fChain->SetBranchAddress("Electron_cutBased_Fall17_V1", &Electron_cutBased_Fall17_V1, &b_Electron_cutBased_Fall17_V1);
   fChain->SetBranchAddress("Electron_jetIdx", &Electron_jetIdx, &b_Electron_jetIdx);
   fChain->SetBranchAddress("Electron_pdgId", &Electron_pdgId, &b_Electron_pdgId);
   fChain->SetBranchAddress("Electron_photonIdx", &Electron_photonIdx, &b_Electron_photonIdx);
   fChain->SetBranchAddress("Electron_tightCharge", &Electron_tightCharge, &b_Electron_tightCharge);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmap", &Electron_vidNestedWPBitmap, &b_Electron_vidNestedWPBitmap);
   fChain->SetBranchAddress("Electron_convVeto", &Electron_convVeto, &b_Electron_convVeto);
   fChain->SetBranchAddress("Electron_cutBased_HEEP", &Electron_cutBased_HEEP, &b_Electron_cutBased_HEEP);
   fChain->SetBranchAddress("Electron_isPFcand", &Electron_isPFcand, &b_Electron_isPFcand);
   fChain->SetBranchAddress("Electron_lostHits", &Electron_lostHits, &b_Electron_lostHits);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WP80", &Electron_mvaFall17V1noIso_WP80, &b_Electron_mvaFall17V1noIso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WP90", &Electron_mvaFall17V1noIso_WP90, &b_Electron_mvaFall17V1noIso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WPL", &Electron_mvaFall17V1noIso_WPL, &b_Electron_mvaFall17V1noIso_WPL);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WP80", &Electron_mvaFall17V2noIso_WP80, &b_Electron_mvaFall17V2noIso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WP90", &Electron_mvaFall17V2noIso_WP90, &b_Electron_mvaFall17V2noIso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WPL", &Electron_mvaFall17V2noIso_WPL, &b_Electron_mvaFall17V2noIso_WPL);
   fChain->SetBranchAddress("Electron_genPartIdx", &Electron_genPartIdx, &b_Electron_genPartIdx);
   fChain->SetBranchAddress("Electron_genPartFlav", &Electron_genPartFlav, &b_Electron_genPartFlav);
   fChain->SetBranchAddress("Muon_dxy", &Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dxyErr", &Muon_dxyErr, &b_Muon_dxyErr);
   fChain->SetBranchAddress("Muon_dz", &Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_dzErr", &Muon_dzErr, &b_Muon_dzErr);
   fChain->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_ip3d", &Muon_ip3d, &b_Muon_ip3d);
   fChain->SetBranchAddress("Muon_jetRelIso", &Muon_jetRelIso, &b_Muon_jetRelIso);
   fChain->SetBranchAddress("Muon_mass", &Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_miniPFRelIso_all", &Muon_miniPFRelIso_all, &b_Muon_miniPFRelIso_all);
   fChain->SetBranchAddress("Muon_miniPFRelIso_chg", &Muon_miniPFRelIso_chg, &b_Muon_miniPFRelIso_chg);
   fChain->SetBranchAddress("Muon_pfRelIso03_all", &Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
   fChain->SetBranchAddress("Muon_pfRelIso03_chg", &Muon_pfRelIso03_chg, &b_Muon_pfRelIso03_chg);
   fChain->SetBranchAddress("Muon_pfRelIso04_all", &Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("Muon_phi", &Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt", &Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_ptErr", &Muon_ptErr, &b_Muon_ptErr);
   fChain->SetBranchAddress("Muon_segmentComp", &Muon_segmentComp, &b_Muon_segmentComp);
   fChain->SetBranchAddress("Muon_sip3d", &Muon_sip3d, &b_Muon_sip3d);
   fChain->SetBranchAddress("Muon_mvaTTH", &Muon_mvaTTH, &b_Muon_mvaTTH);
   fChain->SetBranchAddress("Muon_charge", &Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_jetIdx", &Muon_jetIdx, &b_Muon_jetIdx);
   fChain->SetBranchAddress("Muon_nStations", &Muon_nStations, &b_Muon_nStations);
   fChain->SetBranchAddress("Muon_nTrackerLayers", &Muon_nTrackerLayers, &b_Muon_nTrackerLayers);
   fChain->SetBranchAddress("Muon_pdgId", &Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_tightCharge", &Muon_tightCharge, &b_Muon_tightCharge);
   fChain->SetBranchAddress("Muon_highPtId", &Muon_highPtId, &b_Muon_highPtId);
   fChain->SetBranchAddress("Muon_inTimeMuon", &Muon_inTimeMuon, &b_Muon_inTimeMuon);
   fChain->SetBranchAddress("Muon_isGlobal", &Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_isPFcand", &Muon_isPFcand, &b_Muon_isPFcand);
   fChain->SetBranchAddress("Muon_isTracker", &Muon_isTracker, &b_Muon_isTracker);
   fChain->SetBranchAddress("Muon_mediumId", &Muon_mediumId, &b_Muon_mediumId);
   fChain->SetBranchAddress("Muon_mediumPromptId", &Muon_mediumPromptId, &b_Muon_mediumPromptId);
   fChain->SetBranchAddress("Muon_miniIsoId", &Muon_miniIsoId, &b_Muon_miniIsoId);
   fChain->SetBranchAddress("Muon_multiIsoId", &Muon_multiIsoId, &b_Muon_multiIsoId);
   fChain->SetBranchAddress("Muon_mvaId", &Muon_mvaId, &b_Muon_mvaId);
   fChain->SetBranchAddress("Muon_pfIsoId", &Muon_pfIsoId, &b_Muon_pfIsoId);
   fChain->SetBranchAddress("Muon_softId", &Muon_softId, &b_Muon_softId);
   fChain->SetBranchAddress("Muon_softMvaId", &Muon_softMvaId, &b_Muon_softMvaId);
   fChain->SetBranchAddress("Muon_tightId", &Muon_tightId, &b_Muon_tightId);
   fChain->SetBranchAddress("Muon_tkIsoId", &Muon_tkIsoId, &b_Muon_tkIsoId);
   fChain->SetBranchAddress("Muon_triggerIdLoose", &Muon_triggerIdLoose, &b_Muon_triggerIdLoose);

   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("weight",1);
   fChain->SetBranchStatus("MET_pt",1);
   fChain->SetBranchStatus("nGenPart",1);
   fChain->SetBranchStatus("GenPart_eta",1);
   fChain->SetBranchStatus("GenPart_mass",1);
   fChain->SetBranchStatus("GenPart_phi",1);
   fChain->SetBranchStatus("GenPart_pt",1);
   fChain->SetBranchStatus("GenPart_genPartIdxMother",1);
   fChain->SetBranchStatus("GenPart_pdgId",1);
   fChain->SetBranchStatus("GenPart_status",1);
   fChain->SetBranchStatus("Nele",1);
   fChain->SetBranchStatus("PT_ele",1);
   fChain->SetBranchStatus("Eta_ele",1);
   fChain->SetBranchStatus("Phi_ele",1);
   fChain->SetBranchStatus("M_ele",1);
   fChain->SetBranchStatus("Electron_pt",1);
   fChain->SetBranchStatus("ID_ele",1);
   fChain->SetBranchStatus("MiniIso_ele",1);
   fChain->SetBranchStatus("FO_baseline_ele",1);
   fChain->SetBranchStatus("baseline_loose_ele",1);
   fChain->SetBranchStatus("baseline_tight_ele",1);
   fChain->SetBranchStatus("Electron_mvaFall17V1noIso_WP80",1);
   fChain->SetBranchStatus("Electron_mvaFall17V1noIso_WP90",1);
   fChain->SetBranchStatus("Electron_mvaFall17V1noIso_WPL",1);
   fChain->SetBranchStatus("Electron_mvaFall17V2noIso_WP80",1);
   fChain->SetBranchStatus("Electron_mvaFall17V2noIso_WP90",1);
   fChain->SetBranchStatus("Electron_mvaFall17V2noIso_WPL",1);
   fChain->SetBranchStatus("Electron_lostHits",1);
   fChain->SetBranchStatus("Electron_ip3d",1);
   fChain->SetBranchStatus("Electron_sip3d",1);
   fChain->SetBranchStatus("Electron_convVeto",1);
   fChain->SetBranchStatus("Electron_dxy",1);
   fChain->SetBranchStatus("Electron_dz",1);
   fChain->SetBranchStatus("genNele",1);
   fChain->SetBranchStatus("genPDGID_ele",1);
   fChain->SetBranchStatus("genMomPDGID_ele",1);
   fChain->SetBranchStatus("genPT_ele",1);
   fChain->SetBranchStatus("genEta_ele",1);
   fChain->SetBranchStatus("genPhi_ele",1);
   fChain->SetBranchStatus("genM_ele",1);
   fChain->SetBranchStatus("Nmu",1);
   fChain->SetBranchStatus("PT_mu",1);
   fChain->SetBranchStatus("Eta_mu",1);
   fChain->SetBranchStatus("Phi_mu",1);
   fChain->SetBranchStatus("M_mu",1);
   fChain->SetBranchStatus("Muon_pt",1);
   fChain->SetBranchStatus("ID_mu",1);
   fChain->SetBranchStatus("MiniIso_mu",1);
   fChain->SetBranchStatus("baseline_loose_mu",1);
   fChain->SetBranchStatus("baseline_tight_mu",1);
   fChain->SetBranchStatus("Muon_mediumId",1);
   fChain->SetBranchStatus("Muon_mediumPromptId",1);
   fChain->SetBranchStatus("Muon_tightId",1);
   fChain->SetBranchStatus("Muon_softId",1);
   fChain->SetBranchStatus("Muon_softMvaId",1);
   fChain->SetBranchStatus("Muon_ip3d",1);
   fChain->SetBranchStatus("Muon_sip3d",1);
   fChain->SetBranchStatus("Muon_dxy",1);
   fChain->SetBranchStatus("Muon_dz",1);
   fChain->SetBranchStatus("genNmu",1);
   fChain->SetBranchStatus("genPDGID_mu",1);
   fChain->SetBranchStatus("genMomPDGID_mu",1);
   fChain->SetBranchStatus("genPT_mu",1);
   fChain->SetBranchStatus("genEta_mu",1);
   fChain->SetBranchStatus("genPhi_mu",1);
   fChain->SetBranchStatus("genM_mu",1);
}

