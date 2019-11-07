#include "ElectronNtuple.hh"
#include "ParticleList.hh"
#include "TInterpreter.h"

//#include "SUSYNANOBase.hh"

ElectronNtuple::ElectronNtuple(TTree* tree)
  : NtupleBase<SUSYNANOBase>(tree){}

ElectronNtuple::~ElectronNtuple(){}

TTree* ElectronNtuple::InitOutputTree(const string& sample){

  gInterpreter->GenerateDictionary("vector<vector<int>>", "vector");
  
  TTree* tree = (TTree*) new TTree(sample.c_str(), sample.c_str());

  tree->Branch("weight", &m_weight);
  
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

  tree->Branch("genNele", &m_genNele);

  tree->Branch("genPT_ele",  &m_genPT_ele);
  tree->Branch("genEta_ele", &m_genEta_ele);
  tree->Branch("genPhi_ele", &m_genPhi_ele);
  tree->Branch("genM_ele",   &m_genM_ele);
  tree->Branch("genCharge_ele",  &m_genCharge_ele);
  tree->Branch("genPDGID_ele",   &m_genPDGID_ele);
  tree->Branch("genMomPDGID_ele",   &m_genMomPDGID_ele);
  tree->Branch("genIndex_ele",   &m_genIndex_ele);

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
  tree->Branch("Electron_mvaSpring16GP", &m_Electron_mvaSpring16GP);
  tree->Branch("Electron_mvaSpring16HZZ", &m_Electron_mvaSpring16HZZ);
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
  tree->Branch("Electron_cutBased_HLTPreSel", &m_Electron_cutBased_HLTPreSel);
  tree->Branch("Electron_cutBased_Spring15", &m_Electron_cutBased_Spring15);
  tree->Branch("Electron_cutBased_Sum16", &m_Electron_cutBased_Sum16);
  tree->Branch("Electron_jetIdx", &m_Electron_jetIdx);
  tree->Branch("Electron_pdgId", &m_Electron_pdgId);
  tree->Branch("Electron_photonIdx", &m_Electron_photonIdx);
  tree->Branch("Electron_tightCharge", &m_Electron_tightCharge);
  tree->Branch("Electron_vidNestedWPBitmap", &m_Electron_vidNestedWPBitmap);
  tree->Branch("Electron_vidNestedWPBitmapSpring15", &m_Electron_vidNestedWPBitmapSpring15);
  tree->Branch("Electron_vidNestedWPBitmapSum16", &m_Electron_vidNestedWPBitmapSum16);
  tree->Branch("Electron_convVeto", &m_Electron_convVeto);
  tree->Branch("Electron_cutBased_HEEP", &m_Electron_cutBased_HEEP);
  tree->Branch("Electron_isPFcand", &m_Electron_isPFcand);
  tree->Branch("Electron_lostHits", &m_Electron_lostHits);
  tree->Branch("Electron_mvaFall17V1Iso_WP80", &m_Electron_mvaFall17V1Iso_WP80);
  tree->Branch("Electron_mvaFall17V1Iso_WP90", &m_Electron_mvaFall17V1Iso_WP90);
  tree->Branch("Electron_mvaFall17V1Iso_WPL", &m_Electron_mvaFall17V1Iso_WPL);
  tree->Branch("Electron_mvaFall17V1noIso_WP80", &m_Electron_mvaFall17V1noIso_WP80);
  tree->Branch("Electron_mvaFall17V1noIso_WP90", &m_Electron_mvaFall17V1noIso_WP90);
  tree->Branch("Electron_mvaFall17V1noIso_WPL", &m_Electron_mvaFall17V1noIso_WPL);
  tree->Branch("Electron_mvaFall17V2Iso_WP80", &m_Electron_mvaFall17V2Iso_WP80);
  tree->Branch("Electron_mvaFall17V2Iso_WP90", &m_Electron_mvaFall17V2Iso_WP90);
  tree->Branch("Electron_mvaFall17V2Iso_WPL", &m_Electron_mvaFall17V2Iso_WPL);
  tree->Branch("Electron_mvaFall17V2noIso_WP80", &m_Electron_mvaFall17V2noIso_WP80);
  tree->Branch("Electron_mvaFall17V2noIso_WP90", &m_Electron_mvaFall17V2noIso_WP90);
  tree->Branch("Electron_mvaFall17V2noIso_WPL", &m_Electron_mvaFall17V2noIso_WPL);
  tree->Branch("Electron_mvaSpring16GP_WP80", &m_Electron_mvaSpring16GP_WP80);
  tree->Branch("Electron_mvaSpring16GP_WP90", &m_Electron_mvaSpring16GP_WP90);
  tree->Branch("Electron_mvaSpring16HZZ_WPL", &m_Electron_mvaSpring16HZZ_WPL);

  tree->Branch("Electron_genPartIdx", &m_Electron_genPartIdx);
  tree->Branch("Electron_genPartFlav", &m_Electron_genPartFlav);

  tree->Branch("Electron_cleanmask", &m_Electron_cleanmask);

  return tree;
}

void ElectronNtuple::ClearVariables(){

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
  m_Electron_mvaSpring16GP.clear();
  m_Electron_mvaSpring16HZZ.clear();
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
  m_Electron_cutBased_HLTPreSel.clear();
  m_Electron_cutBased_Spring15.clear();
  m_Electron_cutBased_Sum16.clear();
  m_Electron_jetIdx.clear();
  m_Electron_pdgId.clear();
  m_Electron_photonIdx.clear();
  m_Electron_tightCharge.clear();
  m_Electron_vidNestedWPBitmap.clear();
  m_Electron_vidNestedWPBitmapSpring15.clear();
  m_Electron_vidNestedWPBitmapSum16.clear();
  m_Electron_convVeto.clear();
  m_Electron_cutBased_HEEP.clear();
  m_Electron_isPFcand.clear();
  m_Electron_lostHits.clear();
  m_Electron_mvaFall17V1Iso_WP80.clear();
  m_Electron_mvaFall17V1Iso_WP90.clear();
  m_Electron_mvaFall17V1Iso_WPL.clear();
  m_Electron_mvaFall17V1noIso_WP80.clear();
  m_Electron_mvaFall17V1noIso_WP90.clear();
  m_Electron_mvaFall17V1noIso_WPL.clear();
  m_Electron_mvaFall17V2Iso_WP80.clear();
  m_Electron_mvaFall17V2Iso_WP90.clear();
  m_Electron_mvaFall17V2Iso_WPL.clear();
  m_Electron_mvaFall17V2noIso_WP80.clear();
  m_Electron_mvaFall17V2noIso_WP90.clear();
  m_Electron_mvaFall17V2noIso_WPL.clear();
  m_Electron_mvaSpring16GP_WP80.clear();
  m_Electron_mvaSpring16GP_WP90.clear();
  m_Electron_mvaSpring16HZZ_WPL.clear();

  m_Electron_genPartIdx.clear();
  m_Electron_genPartFlav.clear();

  m_Electron_cleanmask.clear();

}

void ElectronNtuple::FillOutputTree(TTree* tree){

  ClearVariables();
    
  ParticleList Electrons = AnalysisBase<SUSYNANOBase>::GetElectrons();
  Electrons = Electrons.PtEtaCut(3.5);
  
  Electrons.SortByPt();

  m_Nele = Electrons.size();

  ParticleList GenElectrons = AnalysisBase<SUSYNANOBase>::GetGenElectrons();
  GenElectrons.SortByPt();

  m_genNele = GenElectrons.size();
  
  // Fill reconstructed lepton branches
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
  vector<int> genmatch;
  for(int i = 0; i < m_genNele; i++)
    genmatch.push_back(-1);
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

    //Fill nano SUSY branches
    for(int j =0; j < nElectron; j++){
      TLorentzVector nano_e;
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
	m_Electron_mvaSpring16GP.push_back(Electron_mvaSpring16GP[j]);
	m_Electron_mvaSpring16HZZ.push_back(Electron_mvaSpring16HZZ[j]);
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
	m_Electron_cutBased_HLTPreSel.push_back(Electron_cutBased_HLTPreSel[j]);
	m_Electron_cutBased_Spring15.push_back(Electron_cutBased_Spring15[j]);
	m_Electron_cutBased_Sum16.push_back(Electron_cutBased_Sum16[j]);
	m_Electron_jetIdx.push_back(Electron_jetIdx[j]);
	m_Electron_pdgId.push_back(Electron_pdgId[j]);
	m_Electron_photonIdx.push_back(Electron_photonIdx[j]);
	m_Electron_tightCharge.push_back(Electron_tightCharge[j]);
	m_Electron_vidNestedWPBitmap.push_back(Electron_vidNestedWPBitmap[j]);
	m_Electron_vidNestedWPBitmapSpring15.push_back(Electron_vidNestedWPBitmapSpring15[j]);
	m_Electron_vidNestedWPBitmapSum16.push_back(Electron_vidNestedWPBitmapSum16[j]);
	m_Electron_convVeto.push_back(Electron_convVeto[j]);
	m_Electron_cutBased_HEEP.push_back(Electron_cutBased_HEEP[j]);
	m_Electron_isPFcand.push_back(Electron_isPFcand[j]);
	m_Electron_lostHits.push_back(Electron_lostHits[j]);
	m_Electron_mvaFall17V1Iso_WP80.push_back(Electron_mvaFall17V1Iso_WP80[j]);
	m_Electron_mvaFall17V1Iso_WP90.push_back(Electron_mvaFall17V1Iso_WP90[j]);
	m_Electron_mvaFall17V1Iso_WPL.push_back(Electron_mvaFall17V1Iso_WPL[j]);
	m_Electron_mvaFall17V1noIso_WP80.push_back(Electron_mvaFall17V1noIso_WP80[j]);
	m_Electron_mvaFall17V1noIso_WP90.push_back(Electron_mvaFall17V1noIso_WP90[j]);
	m_Electron_mvaFall17V1noIso_WPL.push_back(Electron_mvaFall17V1noIso_WPL[j]);
	m_Electron_mvaFall17V2Iso_WP80.push_back(Electron_mvaFall17V2Iso_WP80[j]);
	m_Electron_mvaFall17V2Iso_WP90.push_back(Electron_mvaFall17V2Iso_WP90[j]);
	m_Electron_mvaFall17V2Iso_WPL.push_back(Electron_mvaFall17V2Iso_WPL[j]);
	m_Electron_mvaFall17V2noIso_WP80.push_back(Electron_mvaFall17V2noIso_WP80[j]);
	m_Electron_mvaFall17V2noIso_WP90.push_back(Electron_mvaFall17V2noIso_WP90[j]);
	m_Electron_mvaFall17V2noIso_WPL.push_back(Electron_mvaFall17V2noIso_WPL[j]);
	m_Electron_mvaSpring16GP_WP80.push_back(Electron_mvaSpring16GP_WP80[j]);
	m_Electron_mvaSpring16GP_WP90.push_back(Electron_mvaSpring16GP_WP90[j]);
	m_Electron_mvaSpring16HZZ_WPL.push_back(Electron_mvaSpring16HZZ_WPL[j]);
	m_Electron_genPartIdx.push_back(Electron_genPartIdx[j]);
	m_Electron_genPartFlav.push_back(Electron_genPartFlav[j]);
	m_Electron_cleanmask.push_back(Electron_cleanmask[j]);
	break;
      }
      
    }

    int index = -1;
    for(int g = 0; g < m_genNele; g++)
      if(Electrons[r].DeltaR(GenElectrons[g]) < 0.02){
	index = g;
	genmatch[g] = r;
	break;
      }
    m_Index_ele.push_back(index);
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
    m_genIndex_ele.push_back(genmatch[g]);
  }

  // Fill output tree
  if(tree)
    tree->Fill();
  
}

