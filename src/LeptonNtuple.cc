#include "LeptonNtuple.hh"
#include "ParticleList.hh"
#include "TInterpreter.h"

//#include "SUSYNANOBase.hh"

LeptonNtuple::LeptonNtuple(TTree* tree)
  : NtupleBase<SUSYNANOBase>(tree){}

LeptonNtuple::~LeptonNtuple(){}

TTree* LeptonNtuple::InitOutputTree(const string& sample){

  gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
  
  TTree* tree = (TTree*) new TTree(sample.c_str(), sample.c_str());

  tree->Branch("weight", &m_weight);
  
  //electrons
  tree->Branch("Nele", &m_Nele);
  tree->Branch("PT_ele",  &m_PT_ele);
  tree->Branch("Eta_ele", &m_Eta_ele);
  tree->Branch("Phi_ele", &m_Phi_ele);
  tree->Branch("M_ele",   &m_M_ele);
  tree->Branch("Charge_ele",  &m_Charge_ele);
  tree->Branch("PDGID_ele",   &m_PDGID_ele);
  tree->Branch("RelIso_ele",  &m_RelIso_ele);
  tree->Branch("MiniIso_ele", &m_MiniIso_ele);
  tree->Branch("ID_ele",      &m_ID_ele);
  tree->Branch("Index_ele",   &m_Index_ele);
  tree->Branch("FO_baseline_ele", &m_baseline_ele);
  tree->Branch("baseline_loose_ele", &m_baseline_loose_ele);
  tree->Branch("baseline_tight_ele", &m_baseline_tight_ele);

  //gen electrons
  tree->Branch("genNele", &m_genNele);
  tree->Branch("genPT_ele",  &m_genPT_ele);
  tree->Branch("genEta_ele", &m_genEta_ele);
  tree->Branch("genPhi_ele", &m_genPhi_ele);
  tree->Branch("genM_ele",   &m_genM_ele);
  tree->Branch("genCharge_ele",  &m_genCharge_ele);
  tree->Branch("genPDGID_ele",   &m_genPDGID_ele);
  tree->Branch("genMomPDGID_ele",   &m_genMomPDGID_ele);
  tree->Branch("genIndex_ele",   &m_genIndex_ele);

  //muons
  tree->Branch("Nmu", &m_Nmu);
  tree->Branch("PT_mu",  &m_PT_mu);
  tree->Branch("Eta_mu", &m_Eta_mu);
  tree->Branch("Phi_mu", &m_Phi_mu);
  tree->Branch("M_mu",   &m_M_mu);
  tree->Branch("Charge_mu",  &m_Charge_mu);
  tree->Branch("PDGID_mu",   &m_PDGID_mu);
  tree->Branch("RelIso_mu",  &m_RelIso_mu);
  tree->Branch("MiniIso_mu", &m_MiniIso_mu);
  tree->Branch("ID_mu",      &m_ID_mu);
  tree->Branch("Index_mu",   &m_Index_mu);
  tree->Branch("baseline_loose_mu", &m_baseline_loose_mu);
  tree->Branch("baseline_tight_mu", &m_baseline_tight_mu);

  //gen muons
  tree->Branch("genNmu", &m_genNmu);
  tree->Branch("genPT_mu",  &m_genPT_mu);
  tree->Branch("genEta_mu", &m_genEta_mu);
  tree->Branch("genPhi_mu", &m_genPhi_mu);
  tree->Branch("genM_mu",   &m_genM_mu);
  tree->Branch("genCharge_mu",  &m_genCharge_mu);
  tree->Branch("genPDGID_mu",   &m_genPDGID_mu);
  tree->Branch("genMomPDGID_mu",   &m_genMomPDGID_mu);
  tree->Branch("genIndex_mu",   &m_genIndex_mu);

  //NanoSUSY ntuple 
  tree->Branch("Electron_deltaEtaSC",&m_Electron_deltaEtaSC);
  tree->Branch("Electron_dr03EcalRecHitSumEt", &m_Electron_dr03EcalRecHitSumEt);
  tree->Branch("Electron_dr03HcalDepth1TowerSumEt", &m_Electron_dr03HcalDepth1TowerSumEt);
  tree->Branch("Electron_dr03TkSumPt", &m_Electron_dr03TkSumPt);
  tree->Branch("Electron_dr03TkSumPtHEEP", &m_Electron_dr03TkSumPtHEEP);
  tree->Branch("Electron_dxy", &m_Electron_dxy);
  tree->Branch("Electron_dxyErr", &m_Electron_dxyErr);
  tree->Branch("Electron_dz", &m_Electron_dz);
  tree->Branch("Electron_dzErr", &m_Electron_dzErr);
  tree->Branch("Electron_eCorr", &m_Electron_eCorr);
  tree->Branch("Electron_eInvMinusPInv", &m_Electron_eInvMinusPInv);
  tree->Branch("Electron_energyErr", &m_Electron_energyErr);
  tree->Branch("Electron_eta", &m_Electron_eta);
  tree->Branch("Electron_hoe", &m_Electron_hoe);
  tree->Branch("Electron_ip3d", &m_Electron_ip3d);
  tree->Branch("Electron_jetRelIso", &m_Electron_jetRelIso);
  tree->Branch("Electron_mass", &m_Electron_mass);
  tree->Branch("Electron_miniPFRelIso_all", &m_Electron_miniPFRelIso_all);
  tree->Branch("Electron_miniPFRelIso_chg", &m_Electron_miniPFRelIso_chg);
  tree->Branch("Electron_mvaFall17V1Iso", &m_Electron_mvaFall17V1Iso);
  tree->Branch("Electron_mvaFall17V1noIso", &m_Electron_mvaFall17V1noIso);
  tree->Branch("Electron_mvaFall17V2Iso", &m_Electron_mvaFall17V2Iso);
  tree->Branch("Electron_mvaFall17V2noIso", &m_Electron_mvaFall17V2noIso);
  // tree->Branch("Electron_mvaSpring16GP", &m_Electron_mvaSpring16GP);
  // tree->Branch("Electron_mvaSpring16HZZ", &m_Electron_mvaSpring16HZZ);
  tree->Branch("Electron_pfRelIso03_all", &m_Electron_pfRelIso03_all);
  tree->Branch("Electron_pfRelIso03_chg", &m_Electron_pfRelIso03_chg);
  tree->Branch("Electron_phi", &m_Electron_phi);
  tree->Branch("Electron_pt", &m_Electron_pt);
  tree->Branch("Electron_r9", &m_Electron_r9);
  tree->Branch("Electron_sieie", &m_Electron_sieie);
  tree->Branch("Electron_sip3d", &m_Electron_sip3d);
  tree->Branch("Electron_mvaTTH", &m_Electron_mvaTTH);
  tree->Branch("Electron_charge", &m_Electron_charge);
  tree->Branch("Electron_cutBased", &m_Electron_cutBased);
  tree->Branch("Electron_cutBased_Fall17_V1", &m_Electron_cutBased_Fall17_V1);
  //tree->Branch("Electron_cutBased_HLTPreSel", &m_Electron_cutBased_HLTPreSel);
  //tree->Branch("Electron_cutBased_Spring15", &m_Electron_cutBased_Spring15);
  //tree->Branch("Electron_cutBased_Sum16", &m_Electron_cutBased_Sum16);
  tree->Branch("Electron_jetIdx", &m_Electron_jetIdx);
  tree->Branch("Electron_pdgId", &m_Electron_pdgId);
  tree->Branch("Electron_photonIdx", &m_Electron_photonIdx);
  tree->Branch("Electron_tightCharge", &m_Electron_tightCharge);
  tree->Branch("Electron_vidNestedWPBitmap", &m_Electron_vidNestedWPBitmap);
  //tree->Branch("Electron_vidNestedWPBitmapSpring15", &m_Electron_vidNestedWPBitmapSpring15);
  //tree->Branch("Electron_vidNestedWPBitmapSum16", &m_Electron_vidNestedWPBitmapSum16);
  tree->Branch("Electron_convVeto", &m_Electron_convVeto);
  tree->Branch("Electron_cutBased_HEEP", &m_Electron_cutBased_HEEP);
  tree->Branch("Electron_isPFcand", &m_Electron_isPFcand);
  tree->Branch("Electron_lostHits", &m_Electron_lostHits);
  tree->Branch("Electron_mvaFall17V1noIso_WP80", &m_Electron_mvaFall17V1noIso_WP80);
  tree->Branch("Electron_mvaFall17V1noIso_WP90", &m_Electron_mvaFall17V1noIso_WP90);
  tree->Branch("Electron_mvaFall17V1noIso_WPL", &m_Electron_mvaFall17V1noIso_WPL);
  tree->Branch("Electron_mvaFall17V2noIso_WP80", &m_Electron_mvaFall17V2noIso_WP80);
  tree->Branch("Electron_mvaFall17V2noIso_WP90", &m_Electron_mvaFall17V2noIso_WP90);
  tree->Branch("Electron_mvaFall17V2noIso_WPL", &m_Electron_mvaFall17V2noIso_WPL);
  //tree->Branch("Electron_mvaSpring16GP_WP80", &m_Electron_mvaSpring16GP_WP80);
  //tree->Branch("Electron_mvaSpring16GP_WP90", &m_Electron_mvaSpring16GP_WP90);
  //tree->Branch("Electron_mvaSpring16HZZ_WPL", &m_Electron_mvaSpring16HZZ_WPL);
  tree->Branch("Electron_genPartIdx", &m_Electron_genPartIdx);
  tree->Branch("Electron_genPartFlav", &m_Electron_genPartFlav);

  tree->Branch("Muon_dxy", &m_Muon_dxy);
  tree->Branch("Muon_dxyErr", &m_Muon_dxyErr);
  tree->Branch("Muon_dz", &m_Muon_dz);
  tree->Branch("Muon_dzErr", &m_Muon_dzErr);
  tree->Branch("Muon_eta", &m_Muon_eta);
  tree->Branch("Muon_ip3d", &m_Muon_ip3d);
  tree->Branch("Muon_jetRelIso", &m_Muon_jetRelIso);
  tree->Branch("Muon_mass", &m_Muon_mass);
  tree->Branch("Muon_miniPFRelIso_all", &m_Muon_miniPFRelIso_all);
  tree->Branch("Muon_miniPFRelIso_chg", &m_Muon_miniPFRelIso_chg);
  tree->Branch("Muon_pfRelIso03_all", &m_Muon_pfRelIso03_all);
  tree->Branch("Muon_pfRelIso03_chg", &m_Muon_pfRelIso03_chg);
  tree->Branch("Muon_pfRelIso04_all", &m_Muon_pfRelIso04_all);
  tree->Branch("Muon_phi", &m_Muon_phi);
  tree->Branch("Muon_pt", &m_Muon_pt);
  tree->Branch("Muon_ptErr", &m_Muon_ptErr);
  tree->Branch("Muon_segmentComp", &m_Muon_segmentComp);
  tree->Branch("Muon_sip3d", &m_Muon_sip3d);
  tree->Branch("Muon_mvaTTH", &m_Muon_mvaTTH);
  tree->Branch("Muon_charge", &m_Muon_charge);
  tree->Branch("Muon_jetIdx", &m_Muon_jetIdx);
  tree->Branch("Muon_nStations", &m_Muon_nStations);
  tree->Branch("Muon_nTrackerLayers", &m_Muon_nTrackerLayers);
  tree->Branch("Muon_pdgId", &m_Muon_pdgId);
  tree->Branch("Muon_tightCharge", &m_Muon_tightCharge);
  tree->Branch("Muon_highPtId", &m_Muon_highPtId);
  tree->Branch("Muon_inTimeMuon", &m_Muon_inTimeMuon);
  tree->Branch("Muon_isGlobal", &m_Muon_isGlobal);
  tree->Branch("Muon_isPFcand", &m_Muon_isPFcand);
  tree->Branch("Muon_isTracker", &m_Muon_isTracker);
  tree->Branch("Muon_mediumId", &m_Muon_mediumId);
  tree->Branch("Muon_mediumPromptId", &m_Muon_mediumPromptId);
  tree->Branch("Muon_miniIsoId", &m_Muon_miniIsoId);
  tree->Branch("Muon_multiIsoId", &m_Muon_multiIsoId);
  tree->Branch("Muon_mvaId", &m_Muon_mvaId);
  tree->Branch("Muon_pfIsoId", &m_Muon_pfIsoId);
  tree->Branch("Muon_softId", &m_Muon_softId);
  tree->Branch("Muon_softMvaId", &m_Muon_softMvaId);
  tree->Branch("Muon_tightId", &m_Muon_tightId);
  tree->Branch("Muon_tkIsoId", &m_Muon_tkIsoId);
  tree->Branch("Muon_triggerIdLoose", &m_Muon_triggerIdLoose);

  return tree;
}

void LeptonNtuple::ClearVariables(){

  m_Electron_deltaEtaSC.clear();
  m_Electron_dr03EcalRecHitSumEt.clear();
  m_Electron_dr03HcalDepth1TowerSumEt.clear();
  m_Electron_dr03TkSumPt.clear();
  m_Electron_dr03TkSumPtHEEP.clear();
  m_Electron_dxy.clear();
  m_Electron_dxyErr.clear();
  m_Electron_dz.clear();
  m_Electron_dzErr.clear();
  m_Electron_eCorr.clear();
  m_Electron_eInvMinusPInv.clear();
  m_Electron_energyErr.clear();
  m_Electron_eta.clear();
  m_Electron_hoe.clear();
  m_Electron_ip3d.clear();
  m_Electron_jetRelIso.clear();
  m_Electron_mass.clear();
  m_Electron_miniPFRelIso_all.clear();
  m_Electron_miniPFRelIso_chg.clear();
  m_Electron_mvaFall17V1Iso.clear();
  m_Electron_mvaFall17V1noIso.clear();
  m_Electron_mvaFall17V2Iso.clear();
  m_Electron_mvaFall17V2noIso.clear();
  // m_Electron_mvaSpring16GP.clear();
  // m_Electron_mvaSpring16HZZ.clear();
  m_Electron_pfRelIso03_all.clear();
  m_Electron_pfRelIso03_chg.clear();
  m_Electron_phi.clear();
  m_Electron_pt.clear();
  m_Electron_r9.clear();
  m_Electron_sieie.clear();
  m_Electron_sip3d.clear();
  m_Electron_mvaTTH.clear();
  m_Electron_charge.clear();
  m_Electron_cutBased.clear();
  m_Electron_cutBased_Fall17_V1.clear();
  // m_Electron_cutBased_HLTPreSel.clear();
  // m_Electron_cutBased_Spring15.clear();
  // m_Electron_cutBased_Sum16.clear();
  m_Electron_jetIdx.clear();
  m_Electron_pdgId.clear();
  m_Electron_photonIdx.clear();
  m_Electron_tightCharge.clear();
  m_Electron_vidNestedWPBitmap.clear();
  // m_Electron_vidNestedWPBitmapSpring15.clear();
  // m_Electron_vidNestedWPBitmapSum16.clear();
  m_Electron_convVeto.clear();
  m_Electron_cutBased_HEEP.clear();
  m_Electron_isPFcand.clear();
  m_Electron_lostHits.clear();
  m_Electron_mvaFall17V1noIso_WP80.clear();
  m_Electron_mvaFall17V1noIso_WP90.clear();
  m_Electron_mvaFall17V1noIso_WPL.clear();
  m_Electron_mvaFall17V2noIso_WP80.clear();
  m_Electron_mvaFall17V2noIso_WP90.clear();
  m_Electron_mvaFall17V2noIso_WPL.clear();
  //m_Electron_mvaSpring16GP_WP80.clear();
  //m_Electron_mvaSpring16GP_WP90.clear();
  //m_Electron_mvaSpring16HZZ_WPL.clear();
  m_Electron_genPartIdx.clear();
  m_Electron_genPartFlav.clear();

  m_Muon_dxy.clear();
  m_Muon_dxyErr.clear();
  m_Muon_dz.clear();
  m_Muon_dzErr.clear();
  m_Muon_eta.clear();
  m_Muon_ip3d.clear();
  m_Muon_jetRelIso.clear();
  m_Muon_mass.clear();
  m_Muon_miniPFRelIso_all.clear();
  m_Muon_miniPFRelIso_chg.clear();
  m_Muon_pfRelIso03_all.clear();
  m_Muon_pfRelIso03_chg.clear();
  m_Muon_pfRelIso04_all.clear();
  m_Muon_phi.clear();
  m_Muon_pt.clear();
  m_Muon_ptErr.clear();
  m_Muon_segmentComp.clear();
  m_Muon_sip3d.clear();
  m_Muon_mvaTTH.clear();
  m_Muon_charge.clear();
  m_Muon_jetIdx.clear();
  m_Muon_nStations.clear();
  m_Muon_nTrackerLayers.clear();
  m_Muon_pdgId.clear();
  m_Muon_tightCharge.clear();
  m_Muon_highPtId.clear();
  m_Muon_inTimeMuon.clear();
  m_Muon_isGlobal.clear();
  m_Muon_isPFcand.clear();
  m_Muon_isTracker.clear();
  m_Muon_mediumId.clear();
  m_Muon_mediumPromptId.clear();
  m_Muon_miniIsoId.clear();
  m_Muon_multiIsoId.clear();
  m_Muon_mvaId.clear();
  m_Muon_pfIsoId.clear();
  m_Muon_softId.clear();
  m_Muon_softMvaId.clear();
  m_Muon_tightId.clear();
  m_Muon_tkIsoId.clear();
  m_Muon_triggerIdLoose.clear();

}

void LeptonNtuple::FillOutputTree(TTree* tree){

  //InitNanoSUSYBranches();
  ClearVariables();
    
  //Electrons
  ParticleList Electrons = AnalysisBase<SUSYNANOBase>::GetElectrons();
  Electrons = Electrons.PtEtaCut(3.5);
  Electrons.SortByPt();
  m_Nele = Electrons.size();

  //Gen Electrons
  ParticleList GenElectrons = AnalysisBase<SUSYNANOBase>::GetGenElectrons();
  GenElectrons.SortByPt();
  m_genNele = GenElectrons.size();

  //Muons
  ParticleList Muons = AnalysisBase<SUSYNANOBase>::GetMuons();
  Muons = Muons.PtEtaCut(3.5);
  Muons.SortByPt();
  m_Nmu = Muons.size();

  //Gen Muons
  ParticleList GenMuons = AnalysisBase<SUSYNANOBase>::GetGenMuons();
  GenMuons.SortByPt();
  m_genNmu = GenMuons.size();

  m_weight = AnalysisBase<SUSYNANOBase>::GetEventWeight();

  // Fill reconstructed electron branches
  m_PT_ele.clear();
  m_Eta_ele.clear();
  m_Phi_ele.clear();
  m_M_ele.clear();
  m_Charge_ele.clear();
  m_PDGID_ele.clear();
  m_RelIso_ele.clear();
  m_MiniIso_ele.clear();
  m_ID_ele.clear();
  m_Index_ele.clear();
  m_baseline_loose_ele.clear();
  m_baseline_tight_ele.clear();
  m_baseline_ele.clear();
  vector<int> genmatch_e;
  TLorentzVector nano_e;
  for(int i = 0; i < m_genNele; i++)
    genmatch_e.push_back(-1);
  for(int r = 0; r < m_Nele; r++){
    m_PT_ele.push_back(Electrons[r].Pt());
    m_Eta_ele.push_back(Electrons[r].Eta());
    m_Phi_ele.push_back(Electrons[r].Phi());
    m_M_ele.push_back(Electrons[r].M());
    m_Charge_ele.push_back(Electrons[r].Charge());
    m_PDGID_ele.push_back(Electrons[r].PDGID());
    m_RelIso_ele.push_back(Electrons[r].RelIso());
    m_MiniIso_ele.push_back(Electrons[r].MiniIso());
    m_ID_ele.push_back(Electrons[r].ParticleID());

    if(Electron_lostHits[r] == 0 && Electron_convVeto[r])
      m_baseline_ele.push_back(true);
    else
      m_baseline_ele.push_back(false);

    if(fabs(Electron_dxy[r]) < 0.05 && fabs(Electron_dz[r]) < 0.1){
      if(Electron_ip3d[r] < 0.0175 && Electron_sip3d[r] < 2.5)
	m_baseline_loose_ele.push_back(true);
      else
	m_baseline_loose_ele.push_back(false);

      if(Electron_ip3d[r] < 0.01 && Electron_sip3d[r] < 2.)
        m_baseline_tight_ele.push_back(true);
      else
	m_baseline_tight_ele.push_back(false);
      }
    else{
      m_baseline_loose_ele.push_back(false);
      m_baseline_tight_ele.push_back(false);
    }

    //Fill nano SUSY branches
    for(int j =0; j < nElectron; j++){
      nano_e.SetPtEtaPhiM(Electron_pt[j],Electron_eta[j],Electron_phi[j],Electron_mass[j]);
      if(Electrons[r].DeltaR(nano_e) < 0.0001){
	m_Electron_deltaEtaSC.push_back(Electron_deltaEtaSC[j]);
	m_Electron_dr03EcalRecHitSumEt.push_back(Electron_dr03EcalRecHitSumEt[j]);
	m_Electron_dr03HcalDepth1TowerSumEt.push_back(Electron_dr03HcalDepth1TowerSumEt[j]);
	m_Electron_dr03TkSumPt.push_back(Electron_dr03TkSumPt[j]);
	m_Electron_dr03TkSumPtHEEP.push_back(Electron_dr03TkSumPtHEEP[j]);
	m_Electron_dxy.push_back(Electron_dxy[j]);
	m_Electron_dxyErr.push_back(Electron_dxyErr[j]);
	m_Electron_dz.push_back(Electron_dz[j]);
	m_Electron_dzErr.push_back(Electron_dzErr[j]);
	m_Electron_eCorr.push_back(Electron_eCorr[j]);
	m_Electron_eInvMinusPInv.push_back(Electron_eInvMinusPInv[j]);
	m_Electron_energyErr.push_back(Electron_energyErr[j]);
	m_Electron_eta.push_back(Electron_eta[j]);
	m_Electron_hoe.push_back(Electron_hoe[j]);
	m_Electron_ip3d.push_back(Electron_ip3d[j]);
	m_Electron_jetRelIso.push_back(Electron_jetRelIso[j]);
	m_Electron_mass.push_back(Electron_mass[j]);
	m_Electron_miniPFRelIso_all.push_back(Electron_miniPFRelIso_all[j]);
	m_Electron_miniPFRelIso_chg.push_back(Electron_miniPFRelIso_chg[j]);
	m_Electron_mvaFall17V1Iso.push_back(Electron_mvaFall17V1Iso[j]);
	m_Electron_mvaFall17V1noIso.push_back(Electron_mvaFall17V1noIso[j]);
	m_Electron_mvaFall17V2Iso.push_back(Electron_mvaFall17V2Iso[j]);
	m_Electron_mvaFall17V2noIso.push_back(Electron_mvaFall17V2noIso[j]);
	// m_Electron_mvaSpring16GP.push_back(Electron_mvaSpring16GP[j]);
	// m_Electron_mvaSpring16HZZ.push_back(Electron_mvaSpring16HZZ[j]);
	m_Electron_pfRelIso03_all.push_back(Electron_pfRelIso03_all[j]);
	m_Electron_pfRelIso03_chg.push_back(Electron_pfRelIso03_chg[j]);
	m_Electron_phi.push_back(Electron_phi[j]);
	m_Electron_pt.push_back(Electron_pt[j]);
	m_Electron_r9.push_back(Electron_r9[j]);
	m_Electron_sieie.push_back(Electron_sieie[j]);
	m_Electron_sip3d.push_back(Electron_sip3d[j]);
	m_Electron_mvaTTH.push_back(Electron_mvaTTH[j]);
	m_Electron_charge.push_back(Electron_charge[j]);
	m_Electron_cutBased.push_back(Electron_cutBased[j]);
	m_Electron_cutBased_Fall17_V1.push_back(Electron_cutBased_Fall17_V1[j]);
	// m_Electron_cutBased_HLTPreSel.push_back(Electron_cutBased_HLTPreSel[j]);
	// m_Electron_cutBased_Spring15.push_back(Electron_cutBased_Spring15[j]);
	// m_Electron_cutBased_Sum16.push_back(Electron_cutBased_Sum16[j]);
	m_Electron_jetIdx.push_back(Electron_jetIdx[j]);
	m_Electron_pdgId.push_back(Electron_pdgId[j]);
	m_Electron_photonIdx.push_back(Electron_photonIdx[j]);
	m_Electron_tightCharge.push_back(Electron_tightCharge[j]);
	m_Electron_vidNestedWPBitmap.push_back(Electron_vidNestedWPBitmap[j]);
	// m_Electron_vidNestedWPBitmapSpring15.push_back(Electron_vidNestedWPBitmapSpring15[j]);
	// m_Electron_vidNestedWPBitmapSum16.push_back(Electron_vidNestedWPBitmapSum16[j]);
	m_Electron_convVeto.push_back(Electron_convVeto[j]);
	m_Electron_cutBased_HEEP.push_back(Electron_cutBased_HEEP[j]);
	m_Electron_isPFcand.push_back(Electron_isPFcand[j]);
	m_Electron_lostHits.push_back(Electron_lostHits[j]);
	m_Electron_mvaFall17V1noIso_WP80.push_back(Electron_mvaFall17V1noIso_WP80[j]);
	m_Electron_mvaFall17V1noIso_WP90.push_back(Electron_mvaFall17V1noIso_WP90[j]);
	m_Electron_mvaFall17V1noIso_WPL.push_back(Electron_mvaFall17V1noIso_WPL[j]);
	m_Electron_mvaFall17V2noIso_WP80.push_back(Electron_mvaFall17V2noIso_WP80[j]);
	m_Electron_mvaFall17V2noIso_WP90.push_back(Electron_mvaFall17V2noIso_WP90[j]);
	m_Electron_mvaFall17V2noIso_WPL.push_back(Electron_mvaFall17V2noIso_WPL[j]);
	// m_Electron_mvaSpring16GP_WP80.push_back(Electron_mvaSpring16GP_WP80[j]);
	// m_Electron_mvaSpring16GP_WP90.push_back(Electron_mvaSpring16GP_WP90[j]);
	// m_Electron_mvaSpring16HZZ_WPL.push_back(Electron_mvaSpring16HZZ_WPL[j]);
	m_Electron_genPartIdx.push_back(Electron_genPartIdx[j]);
	m_Electron_genPartFlav.push_back(Electron_genPartFlav[j]);
	break;
      }
      
    }

    int index = -1;
    for(int g = 0; g < m_genNele; g++)
      if(Electrons[r].DeltaR(GenElectrons[g]) < 0.02){
	index = g;
	genmatch_e[g] = r;
	break;
      }
    m_Index_ele.push_back(index);
  }

  // Fill reconstructed muon branches
  m_PT_mu.clear();
  m_Eta_mu.clear();
  m_Phi_mu.clear();
  m_M_mu.clear();
  m_Charge_mu.clear();
  m_PDGID_mu.clear();
  m_RelIso_mu.clear();
  m_MiniIso_mu.clear();
  m_ID_mu.clear();
  m_Index_mu.clear();
  m_baseline_loose_mu.clear();
  m_baseline_tight_mu.clear();
  vector<int> genmatch_mu;
  TLorentzVector nano_mu;
  for(int i = 0; i < m_genNmu; i++)
    genmatch_mu.push_back(-1);
  for(int r = 0; r < m_Nmu; r++){
    m_PT_mu.push_back(Muons[r].Pt());
    m_Eta_mu.push_back(Muons[r].Eta());
    m_Phi_mu.push_back(Muons[r].Phi());
    m_M_mu.push_back(Muons[r].M());
    m_Charge_mu.push_back(Muons[r].Charge());
    m_PDGID_mu.push_back(Muons[r].PDGID());
    m_RelIso_mu.push_back(Muons[r].RelIso());
    m_MiniIso_mu.push_back(Muons[r].MiniIso());
    m_ID_mu.push_back(Muons[r].ParticleID());

    if(fabs(Muon_dxy[r]) < 0.05 && fabs(Muon_dz[r]) < 0.1){
      if(Muon_ip3d[r] < 0.0175 && Muon_sip3d[r] < 2.5)
        m_baseline_loose_mu.push_back(true);
      else
	m_baseline_loose_mu.push_back(false);

      if(Muon_ip3d[r] < 0.01 && Muon_sip3d[r] < 2.)
        m_baseline_tight_mu.push_back(true);
      else
	m_baseline_tight_mu.push_back(false);
    }
    else{
      m_baseline_loose_mu.push_back(false);
      m_baseline_tight_mu.push_back(false);
    }

    for(int j =0; j < nMuon; j++){
      nano_mu.SetPtEtaPhiM(Muon_pt[j],Muon_eta[j],Muon_phi[j],Muon_mass[j]);
      if(Muons[r].DeltaR(nano_mu) < 0.0001){
	m_Muon_dxy.push_back(Muon_dxy[j]);
	m_Muon_dxyErr.push_back(Muon_dxyErr[j]);
	m_Muon_dz.push_back(Muon_dz[j]);
	m_Muon_dzErr.push_back(Muon_dzErr[j]);
	m_Muon_eta.push_back(Muon_eta[j]);
	m_Muon_ip3d.push_back(Muon_ip3d[j]);
	m_Muon_jetRelIso.push_back(Muon_jetRelIso[j]);
	m_Muon_mass.push_back(Muon_mass[j]);
	m_Muon_miniPFRelIso_all.push_back(Muon_miniPFRelIso_all[j]);
	m_Muon_miniPFRelIso_chg.push_back(Muon_miniPFRelIso_chg[j]);
	m_Muon_pfRelIso03_all.push_back(Muon_pfRelIso03_all[j]);
	m_Muon_pfRelIso03_chg.push_back(Muon_pfRelIso03_chg[j]);
	m_Muon_pfRelIso04_all.push_back(Muon_pfRelIso04_all[j]);
	m_Muon_phi.push_back(Muon_phi[j]);
	m_Muon_pt.push_back(Muon_pt[j]);
	m_Muon_ptErr.push_back(Muon_ptErr[j]);
	m_Muon_segmentComp.push_back(Muon_segmentComp[j]);
	m_Muon_sip3d.push_back(Muon_sip3d[j]);
	m_Muon_mvaTTH.push_back(Muon_mvaTTH[j]);
	m_Muon_charge.push_back(Muon_charge[j]);
	m_Muon_jetIdx.push_back(Muon_jetIdx[j]);
	m_Muon_nStations.push_back(Muon_nStations[j]);
	m_Muon_nTrackerLayers.push_back(Muon_nTrackerLayers[j]);
	m_Muon_pdgId.push_back(Muon_pdgId[j]);
	m_Muon_tightCharge.push_back(Muon_tightCharge[j]);
	m_Muon_highPtId.push_back(Muon_highPtId[j]);
	m_Muon_inTimeMuon.push_back(Muon_inTimeMuon[j]);
	m_Muon_isGlobal.push_back(Muon_isGlobal[j]);
	m_Muon_isPFcand.push_back(Muon_isPFcand[j]);
	m_Muon_isTracker.push_back(Muon_isTracker[j]);
	m_Muon_mediumId.push_back(Muon_mediumId[j]);
	m_Muon_mediumPromptId.push_back(Muon_mediumPromptId[j]);
	m_Muon_miniIsoId.push_back(Muon_miniIsoId[j]);
	m_Muon_multiIsoId.push_back(Muon_multiIsoId[j]);
	m_Muon_mvaId.push_back(Muon_mvaId[j]);
	m_Muon_pfIsoId.push_back(Muon_pfIsoId[j]);
	m_Muon_softId.push_back(Muon_softId[j]);
	m_Muon_softMvaId.push_back(Muon_softMvaId[j]);
	m_Muon_tightId.push_back(Muon_tightId[j]);
	m_Muon_tkIsoId.push_back(Muon_tkIsoId[j]);
	m_Muon_triggerIdLoose.push_back(Muon_triggerIdLoose[j]);
	break;
      }
    }
    int index = -1;
    for(int g = 0; g < m_genNmu; g++)
      if(Muons[r].DeltaR(GenMuons[g]) < 0.02){
        index = g;
        genmatch_mu[g] = r;
        break;
      }
    m_Index_mu.push_back(index);
  }
  
  // Fill gen lepton branches
  m_genPT_ele.clear();
  m_genEta_ele.clear();
  m_genPhi_ele.clear();
  m_genM_ele.clear();
  m_genCharge_ele.clear();
  m_genPDGID_ele.clear();
  m_genMomPDGID_ele.clear();
  m_genIndex_ele.clear();
  for(int g = 0; g < m_genNele; g++){
    m_genPT_ele.push_back(GenElectrons[g].Pt());
    m_genEta_ele.push_back(GenElectrons[g].Eta());
    m_genPhi_ele.push_back(GenElectrons[g].Phi());
    m_genM_ele.push_back(GenElectrons[g].M());
    m_genCharge_ele.push_back(GenElectrons[g].Charge());
    m_genPDGID_ele.push_back(GenElectrons[g].PDGID());
    m_genMomPDGID_ele.push_back(GenElectrons[g].MomPDGID());
    m_genIndex_ele.push_back(genmatch_e[g]);
  }

  m_genPT_mu.clear();
  m_genEta_mu.clear();
  m_genPhi_mu.clear();
  m_genM_mu.clear();
  m_genCharge_mu.clear();
  m_genPDGID_mu.clear();
  m_genMomPDGID_mu.clear();
  m_genIndex_mu.clear();
  for(int g = 0; g < m_genNmu; g++){
    m_genPT_mu.push_back(GenMuons[g].Pt());
    m_genEta_mu.push_back(GenMuons[g].Eta());
    m_genPhi_mu.push_back(GenMuons[g].Phi());
    m_genM_mu.push_back(GenMuons[g].M());
    m_genCharge_mu.push_back(GenMuons[g].Charge());
    m_genPDGID_mu.push_back(GenMuons[g].PDGID());
    m_genMomPDGID_mu.push_back(GenMuons[g].MomPDGID());
    m_genIndex_mu.push_back(genmatch_mu[g]);
  }

  // Fill output tree
  if(tree)
    tree->Fill();
  
}

void LeptonNtuple::InitNanoSUSYBranches(){

  fChain->SetBranchStatus("*", 0);

  fChain->SetBranchStatus("Electron_deltaEtaSC", 1);
  fChain->SetBranchStatus("Electron_dr03EcalRecHitSumEt", 1);
  fChain->SetBranchStatus("Electron_dr03HcalDepth1TowerSumEt", 1);
  fChain->SetBranchStatus("Electron_dr03TkSumPt", 1);
  fChain->SetBranchStatus("Electron_dr03TkSumPtHEEP", 1);
  fChain->SetBranchStatus("Electron_dxy", 1);
  fChain->SetBranchStatus("Electron_dxyErr", 1);
  fChain->SetBranchStatus("Electron_dz", 1);
  fChain->SetBranchStatus("Electron_dzErr", 1);
  fChain->SetBranchStatus("Electron_eCorr", 1);
  fChain->SetBranchStatus("Electron_eInvMinusPInv", 1); 
  fChain->SetBranchStatus("Electron_energyErr", 1);
  fChain->SetBranchStatus("Electron_eta", 1);
  fChain->SetBranchStatus("Electron_hoe", 1);
  fChain->SetBranchStatus("Electron_ip3d", 1);
  fChain->SetBranchStatus("Electron_jetRelIso", 1);
  fChain->SetBranchStatus("Electron_mass", 1);
  fChain->SetBranchStatus("Electron_miniPFRelIso_all", 1);
  fChain->SetBranchStatus("Electron_miniPFRelIso_chg", 1);
  fChain->SetBranchStatus("Electron_mvaFall17V1Iso",  1);
  fChain->SetBranchStatus("Electron_mvaFall17V1noIso", 1);
  fChain->SetBranchStatus("Electron_mvaFall17V2Iso", 1);
  fChain->SetBranchStatus("Electron_mvaFall17V2noIso", 1);
  // fChain->SetBranchStatus("Electron_mvaSpring16GP", 1);
  // fChain->SetBranchStatus("Electron_mvaSpring16HZZ", 1);
  fChain->SetBranchStatus("Electron_pfRelIso03_all", 1);
  fChain->SetBranchStatus("Electron_pfRelIso03_chg", 1);
  fChain->SetBranchStatus("Electron_phi", 1);
  fChain->SetBranchStatus("Electron_pt", 1);
  fChain->SetBranchStatus("Electron_r9", 1);
  fChain->SetBranchStatus("Electron_sieie", 1);
  fChain->SetBranchStatus("Electron_sip3d", 1);
  fChain->SetBranchStatus("Electron_mvaTTH", 1);
  fChain->SetBranchStatus("Electron_charge", 1);
  fChain->SetBranchStatus("Electron_cutBased", 1);
  fChain->SetBranchStatus("Electron_cutBased_Fall17_V1", 1);
  // fChain->SetBranchStatus("Electron_cutBased_HLTPreSel", 1);
  // fChain->SetBranchStatus("Electron_cutBased_Spring15", 1);
  // fChain->SetBranchStatus("Electron_cutBased_Sum16", 1);
  fChain->SetBranchStatus("Electron_jetIdx", 1);
  fChain->SetBranchStatus("Electron_pdgId", 1);
  fChain->SetBranchStatus("Electron_photonIdx", 1);
  fChain->SetBranchStatus("Electron_tightCharge", 1);
  fChain->SetBranchStatus("Electron_vidNestedWPBitmap", 1);
  // fChain->SetBranchStatus("Electron_vidNestedWPBitmapSpring15", 1);
  // fChain->SetBranchStatus("Electron_vidNestedWPBitmapSum16", 1);
  fChain->SetBranchStatus("Electron_convVeto", 1);
  fChain->SetBranchStatus("Electron_cutBased_HEEP", 1);
  fChain->SetBranchStatus("Electron_isPFcand", 1);
  fChain->SetBranchStatus("Electron_lostHits", 1);
  fChain->SetBranchStatus("Electron_mvaFall17V1noIso_WP80", 1);
  fChain->SetBranchStatus("Electron_mvaFall17V1noIso_WP90", 1);
  fChain->SetBranchStatus("Electron_mvaFall17V1noIso_WPL", 1);
  fChain->SetBranchStatus("Electron_mvaFall17V2noIso_WP80", 1);
  fChain->SetBranchStatus("Electron_mvaFall17V2noIso_WP90", 1);
  fChain->SetBranchStatus("Electron_mvaFall17V2noIso_WPL", 1);
  // fChain->SetBranchStatus("Electron_mvaSpring16GP_WP80", 1);
  // fChain->SetBranchStatus("Electron_mvaSpring16GP_WP90", 1);
  // fChain->SetBranchStatus("Electron_mvaSpring16HZZ_WPL", 1);
  fChain->SetBranchStatus("Electron_genPartIdx", 1);
  fChain->SetBranchStatus("Electron_genPartFlav", 1);

  fChain->SetBranchStatus("Muon_dxy", 1);
  fChain->SetBranchStatus("Muon_dxyErr", 1);
  fChain->SetBranchStatus("Muon_dz", 1);
  fChain->SetBranchStatus("Muon_dzErr", 1);
  fChain->SetBranchStatus("Muon_eta", 1);
  fChain->SetBranchStatus("Muon_ip3d", 1);
  fChain->SetBranchStatus("Muon_jetRelIso", 1);
  fChain->SetBranchStatus("Muon_mass", 1);
  fChain->SetBranchStatus("Muon_miniPFRelIso_all", 1);
  fChain->SetBranchStatus("Muon_miniPFRelIso_chg", 1);
  fChain->SetBranchStatus("Muon_pfRelIso03_all", 1);
  fChain->SetBranchStatus("Muon_pfRelIso03_chg", 1);
  fChain->SetBranchStatus("Muon_pfRelIso04_all", 1);
  fChain->SetBranchStatus("Muon_phi", 1);
  fChain->SetBranchStatus("Muon_pt", 1);
  fChain->SetBranchStatus("Muon_ptErr", 1);
  fChain->SetBranchStatus("Muon_segmentComp", 1);
  fChain->SetBranchStatus("Muon_sip3d", 1);
  fChain->SetBranchStatus("Muon_mvaTTH", 1);
  fChain->SetBranchStatus("Muon_charge", 1);
  fChain->SetBranchStatus("Muon_jetIdx", 1);
  fChain->SetBranchStatus("Muon_nStations", 1);
  fChain->SetBranchStatus("Muon_nTrackerLayers", 1);
  fChain->SetBranchStatus("Muon_pdgId", 1);
  fChain->SetBranchStatus("Muon_tightCharge", 1);
  fChain->SetBranchStatus("Muon_highPtId", 1);
  fChain->SetBranchStatus("Muon_inTimeMuon", 1);
  fChain->SetBranchStatus("Muon_isGlobal", 1);
  fChain->SetBranchStatus("Muon_isPFcand", 1);
  fChain->SetBranchStatus("Muon_isTracker", 1);
  fChain->SetBranchStatus("Muon_mediumId", 1);
  fChain->SetBranchStatus("Muon_mediumPromptId", 1);
  fChain->SetBranchStatus("Muon_miniIsoId", 1);
  fChain->SetBranchStatus("Muon_multiIsoId", 1);
  fChain->SetBranchStatus("Muon_mvaId", 1);
  fChain->SetBranchStatus("Muon_pfIsoId", 1);
  fChain->SetBranchStatus("Muon_softId", 1);
  fChain->SetBranchStatus("Muon_softMvaId", 1);
  fChain->SetBranchStatus("Muon_tightId", 1);
  fChain->SetBranchStatus("Muon_tkIsoId", 1);
  fChain->SetBranchStatus("Muon_triggerIdLoose", 1);

}
