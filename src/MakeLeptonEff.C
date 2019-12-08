#include <stdio.h>  
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>

#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLorentzVector.h>
#include <TLatex.h>

#include "KUAnalysis.h"

using namespace std;

void PlotEff(vector<TEfficiency*> vec_eff, vector<TString> vec_names, TString name, TString xlabel, TString ylabel, TString sample, bool isEff, float ymax = 1.05, float ymin = 0.);
void Plot1D(TH1F* h1, TString name, TString sample);
void Plot2D(TH2F* h2, TString name, TString xlabel, TString ylabel, TString sample);
void RocCurve(vector<TEfficiency*> eff, vector<TEfficiency*> fakes, vector<TString> vec_names, int nbins, TString name, TString sample);
void CMSmark(TString plotTitle);

TString RmSpace(TString str);
TGraphAsymmErrors* RocFromEff(TEfficiency* eff, TEfficiency* fakes, int nbins);

int main(int argc, char *argv[]){

  int opt;
  char* filename;
  char* txtfile;
  char* samplename;

  bool test = false;

  //Command line option definitions  
  while((opt = getopt(argc, argv, "i:s:f:th")) != -1){
    switch(opt)  
      {
	//Input filename option [-i]
      case 'i':
	filename = optarg;
	break;
	
	//Input filename option [-s]
      case 's':
        samplename = optarg;
        break;

      case 'f':
	txtfile = optarg;
	break;

      case 't':
	test = true;
	break;

	//input Help option [-h]    
      case 'h':
	printf("Usage: %s [-f] inputList.txt [-s] sampleName\n", argv[0]);
	return 0;
	
	//If unknown option is given, return error
      case '?':
	fprintf(stderr,"Usage: %s [-i] inputfile \n", argv[0]);
	exit(EXIT_FAILURE);
      }
  }

  TLorentzVector gen_ele, ele_cand, gen_mu, mu_cand;
  vector<TLorentzVector> gen_ele_vec, gen_mu_vec;
  vector<int> passIdx_ele, passIdx_mu;
  int idx = 0;

  int nbins = 40;
  float xmin = 0.;
  float xmax = 100.;
  int nbinsY = 40;
  float ymin = 0.;
  float ymax1 = 2.;
  float ymax2 = 100.;
  float dR_cut_ele = 0.01;
  float dR_cut_mu = 0.01;
  float ip3D = 0.0175;
  float sip3D = 8.0;
  float dxy = 0.05;
  float dz = 0.1;

  float weight = -999.;
  float pt = -999.;
  float miniIso = -999.;

  bool isGenZ, isGenW;
  bool isGenZtoEle, isGenWtoEle;
  bool isGenZtoMu, isGenWtoMu;
  bool isGenTau, isGenTauToEle, isGenTauToMu;

  TString ylabel_eff = "Efficiency";
  TString ylabel_fr = "Fakerate";
  TString ele_label = "Electron p_{T} [GeV]";
  TString mu_label = "Muon p_{T} [GeV]";

  TH1F *dR_ele = new TH1F("dR_ele","dR_ele", 100, 0., 0.1);
  TH1F *dR_mu = new TH1F("dR_mu","dR_mu", 100, 0., 0.1);

  TH2F *miniIsoVsPt_sig = new TH2F("miniIsoVsPt_sig","miniIsoVsPt_sig",nbins,xmin,xmax,nbinsY,ymin,ymax1);
  TH2F *miniIsoPtVsPt_sig = new TH2F("miniIsoPtVsPt_sig","miniIsoPtVsPt_sig",nbins,xmin,xmax,nbinsY,ymin,ymax2);
  TH2F *miniIsoVsPt_fake = new TH2F("miniIsoVsPt_fake","miniIsoVsPt_fake",nbins,xmin,xmax,nbinsY,ymin,ymax1);
  TH2F *miniIsoPtVsPt_fake = new TH2F("miniIsoPtVsPt_fake","miniIsoPtVsPt_fake",nbins,xmin,xmax,nbinsY,ymin,ymax2);

  TH2F *miniIsoVsPt_mu_sig = new TH2F("miniIsoVsPt_mu_sig","miniIsoVsPt_mu_sig",nbins,xmin,xmax,nbinsY,ymin,ymax1);
  TH2F *miniIsoPtVsPt_mu_sig = new TH2F("miniIsoPtVsPt_mu_sig","miniIsoPtVsPt_mu_sig",nbins,xmin,xmax,nbinsY,ymin,ymax2);
  TH2F *miniIsoVsPt_mu_fake = new TH2F("miniIsoVsPt_mu_fake","miniIsoVsPt_mu_fake",nbins,xmin,xmax,nbinsY,ymin,ymax1);
  TH2F *miniIsoPtVsPt_mu_fake = new TH2F("miniIsoPtVsPt_mu_fake","miniIsoPtVsPt_mu_fake",nbins,xmin,xmax,nbinsY,ymin,ymax2);

  //Working points electrons (SUSY)
  TH1F *num_vl_ele = new TH1F("num_vl_ele","num_vl_ele",nbins,xmin,xmax);
  TH1F *num_l_ele = new TH1F("num_l_ele","num_l_ele",nbins,xmin,xmax);
  TH1F *num_m_ele = new TH1F("num_m_ele","num_m_ele",nbins,xmin,xmax);
  TH1F *num_t_ele = new TH1F("num_t_ele","num_t_ele",nbins,xmin,xmax);
  
  //Working points for electrons with baseline criteria (SUSY)
  TH1F *num_vl_bl_ele = new TH1F("num_vl_bl_ele","num_vl_bl_ele",nbins,xmin,xmax);
  TH1F *num_l_bl_ele = new TH1F("num_l_bl_ele","num_l_bl_ele",nbins,xmin,xmax);
  TH1F *num_t_bl_ele = new TH1F("num_t_bl_ele","num_t_bl_ele",nbins,xmin,xmax);

  //Working point fakerate (SUSY)
  TH1F *num_bl_fr_ele = new TH1F("num_bl_fr_ele","num_bl_fr_ele",nbins,xmin,xmax);
  TH1F *num_vl_fr_ele = new TH1F("num_vl_fr_ele","num_vl_fr_ele",nbins,xmin,xmax);
  TH1F *num_l_fr_ele = new TH1F("num_l_fr_ele","num_l_fr_ele",nbins,xmin,xmax);
  TH1F *num_t_fr_ele = new TH1F("num_t_fr_ele","num_t_fr_ele",nbins,xmin,xmax);

  //Working points V1 (EGamma)
  TH1F *num_WP80V1 = new TH1F("num_WP80V1","num_WP80V1",nbins,xmin,xmax);
  TH1F *num_WP90V1 = new TH1F("num_WP90V1","num_WP90V1",nbins,xmin,xmax);
  TH1F *num_WPLV1 = new TH1F("num_WPLV1","num_WPLV1",nbins,xmin,xmax);

  //Working points V1 fakerates (EGamma)
  TH1F *num_WP80V1_fr = new TH1F("num_WP80V1_fr","num_WP80V1_fr",nbins,xmin,xmax);
  TH1F *num_WP90V1_fr = new TH1F("num_WP90V1_fr","num_WP90V1_fr",nbins,xmin,xmax);
  TH1F *num_WPLV1_fr = new TH1F("num_WPLV1_fr","num_WPLV1_fr",nbins,xmin,xmax);

  //Working points V2 (EGamma)
  TH1F *num_WP80V2 = new TH1F("num_WP80V2","num_WP80V2",nbins,xmin,xmax);
  TH1F *num_WP90V2 = new TH1F("num_WP90V2","num_WP90V2",nbins,xmin,xmax);
  TH1F *num_WPLV2 = new TH1F("num_WPLV2","num_WPLV2",nbins,xmin,xmax);

  //Working points V2 fakerates (EGamma)
  TH1F *num_WP80V2_fr = new TH1F("num_WP80V2_fr","num_WP80V2_fr",nbins,xmin,xmax);
  TH1F *num_WP90V2_fr = new TH1F("num_WP90V2_fr","num_WP90V2_fr",nbins,xmin,xmax);
  TH1F *num_WPLV2_fr = new TH1F("num_WPLV2_fr","num_WPLV2_fr",nbins,xmin,xmax);

  //Electron candidate selection parameter histograms
  TH1F *num_lostHits = new TH1F("num_lostHits","num_lostHits",nbins,xmin,xmax);
  TH1F *num_ip3D = new TH1F("num_ip3D","num_ip3D",nbins,xmin,xmax);
  TH1F *num_sip3D = new TH1F("num_sip3D","num_sip3D",nbins,xmin,xmax);
  TH1F *num_sip3D_4 = new TH1F("num_sip3D_4","num_sip3D_4",nbins,xmin,xmax);
  TH1F *num_sip3D_8 = new TH1F("num_sip3D_8","num_sip3D_8",nbins,xmin,xmax);
  TH1F *num_convVeto = new TH1F("num_convVeto","num_convVeto",nbins,xmin,xmax);
  TH1F *num_dxy = new TH1F("num_dxy","num_dxy",nbins,xmin,xmax);
  TH1F *num_dz = new TH1F("num_dz","num_dz",nbins,xmin,xmax);

  //Fakerate
  TH1F *num_lostHits_f = new TH1F("num_lostHits_f","num_lostHits_f",nbins,xmin,xmax);
  TH1F *num_ip3D_f = new TH1F("num_ip3D_f","num_ip3D_f",nbins,xmin,xmax);
  TH1F *num_sip3D_f = new TH1F("num_sip3D_f","num_sip3D_f",nbins,xmin,xmax);
  TH1F *num_sip3D_f_4 = new TH1F("num_sip3D_f_4","num_sip3D_f_4",nbins,xmin,xmax);
  TH1F *num_sip3D_f_8 = new TH1F("num_sip3D_f_8","num_sip3D_f_8",nbins,xmin,xmax);
  TH1F *num_convVeto_f = new TH1F("num_convVeto_f","num_convVeto_f",nbins,xmin,xmax);
  TH1F *num_dxy_f = new TH1F("num_dxy_f","num_dxy_f",nbins,xmin,xmax);
  TH1F *num_dz_f = new TH1F("num_dz_f","num_dz_f",nbins,xmin,xmax);

  //Electron Isolation efficiency histograms
  TH1F* num_SUSY_tight_ele_iso0p1  = new TH1F("num_SUSY_tight_ele_iso0p1","num_SUSY_tight_ele_iso0p1",nbins,xmin,xmax);
  TH1F* num_SUSY_tight_ele_iso0p2  = new TH1F("num_SUSY_tight_ele_iso0p2","num_SUSY_tight_ele_iso0p2",nbins,xmin,xmax);
  TH1F* num_SUSY_tight_ele_iso4  = new TH1F("num_SUSY_tight_ele_iso4","num_SUSY_tight_ele_iso4",nbins,xmin,xmax);
  TH1F* num_SUSY_tight_ele_iso6  = new TH1F("num_SUSY_tight_ele_iso6","num_SUSY_tight_ele_iso6",nbins,xmin,xmax);

  //Electron Isolation fakerate histograms
  TH1F* num_SUSY_tight_ele_iso0p1_fr  = new TH1F("num_SUSY_tight_ele_iso0p1_fr","num_SUSY_tight_ele_iso0p1_fr",nbins,xmin,xmax);
  TH1F* num_SUSY_tight_ele_iso0p2_fr  = new TH1F("num_SUSY_tight_ele_iso0p2_fr","num_SUSY_tight_ele_iso0p2_fr",nbins,xmin,xmax);
  TH1F* num_SUSY_tight_ele_iso4_fr  = new TH1F("num_SUSY_tight_ele_iso4_fr","num_SUSY_tight_ele_iso4_fr",nbins,xmin,xmax);
  TH1F* num_SUSY_tight_ele_iso6_fr  = new TH1F("num_SUSY_tight_ele_iso6_fr","num_SUSY_tight_ele_iso6_fr",nbins,xmin,xmax);

  //Muon working point efficiency (Muon POG)
  TH1F *num_bl_mu = new TH1F("num_bl_mu","num_bl_mu",nbins,xmin,xmax);
  TH1F *num_mp_mu = new TH1F("num_mp_mu","num_mp_mu",nbins,xmin,xmax);
  TH1F *num_m_mu = new TH1F("num_m_mu","num_m_mu",nbins,xmin,xmax);
  TH1F *num_t_mu = new TH1F("num_t_mu","num_t_mu",nbins,xmin,xmax);
  TH1F *num_sid_mu = new TH1F("num_sid_mu","num_sid_mu",nbins,xmin,xmax);
  TH1F *num_smva_mu = new TH1F("num_smva_mu","num_smva_mu",nbins,xmin,xmax);

  //Muon working point fakerate (Muon POG)
  TH1F *num_bl_fr_mu = new TH1F("num_bl_fr_mu","num_bl_fr_mu",nbins,xmin,xmax);
  TH1F *num_mp_fr_mu = new TH1F("num_mp_fr_mu","num_mp_fr_mu",nbins,xmin,xmax);
  TH1F *num_m_fr_mu = new TH1F("num_m_fr_mu","num_m_fr_mu",nbins,xmin,xmax);
  TH1F *num_t_fr_mu = new TH1F("num_t_fr_mu","num_t_fr_mu",nbins,xmin,xmax);
  TH1F *num_sid_fr_mu = new TH1F("num_sid_fr_mu","num_sid_fr_mu",nbins,xmin,xmax);
  TH1F *num_smva_fr_mu = new TH1F("num_smva_fr_mu","num_smva_fr_mu",nbins,xmin,xmax);

  //Muon candidate selection parameter efficiency histograms
  TH1F *num_ip3D_mu = new TH1F("num_ip3D_mu","num_ip3D_mu",nbins,xmin,xmax);
  TH1F *num_sip3D_mu = new TH1F("num_sip3D_mu","num_sip3D_mu",nbins,xmin,xmax);
  TH1F *num_sip3D_mu_4 = new TH1F("num_sip3D_mu_4","num_sip3D_mu_4",nbins,xmin,xmax);
  TH1F *num_sip3D_mu_8 = new TH1F("num_sip3D_mu_8","num_sip3D_mu_8",nbins,xmin,xmax);
  TH1F *num_dxy_mu = new TH1F("num_dxy_mu","num_dxy_mu",nbins,xmin,xmax);
  TH1F *num_dz_mu = new TH1F("num_dz_mu","num_dz_mu",nbins,xmin,xmax);

  //Muon candidate selection parameter fakerate histograms
  TH1F *num_ip3D_fr_mu = new TH1F("num_ip3D_fr_mu","num_ip3D_fr_mu",nbins,xmin,xmax);
  TH1F *num_sip3D_fr_mu = new TH1F("num_sip3D_fr_mu","num_sip3D_fr_mu",nbins,xmin,xmax);
  TH1F *num_sip3D_fr_mu_4 = new TH1F("num_sip3D_fr_mu_4","num_sip3D_fr_mu_4",nbins,xmin,xmax);
  TH1F *num_sip3D_fr_mu_8 = new TH1F("num_sip3D_fr_mu_8","num_sip3D_fr_mu_8",nbins,xmin,xmax);
  TH1F *num_dxy_fr_mu = new TH1F("num_dxy_fr_mu","num_dxy_fr_mu",nbins,xmin,xmax);
  TH1F *num_dz_fr_mu = new TH1F("num_dz_fr_mu","num_dz_fr_mu",nbins,xmin,xmax);

  //Denominator true and fake electrons and muons
  TH1F *den_all_ele = new TH1F("den_all_ele","den_all_ele",nbins,xmin,xmax);
  TH1F *den_fake_ele = new TH1F("den_fake_ele","den_fake_ele",nbins,xmin,xmax);
  TH1F *den_all_mu = new TH1F("den_all_mu","den_all_mu",nbins,xmin,xmax);
  TH1F *den_fake_mu = new TH1F("den_fake_mu","den_fake_mu",nbins,xmin,xmax);

  //Muon Isolation efficiency histograms
  TH1F* num_SUSY_medPrompt_mu_iso0p1  = new TH1F("num_SUSY_medPrompt_mu_iso0p1","num_SUSY_medPrompt_mu_iso0p1",nbins,xmin,xmax);
  TH1F* num_SUSY_medPrompt_mu_iso0p2  = new TH1F("num_SUSY_medPrompt_mu_iso0p2","num_SUSY_medPrompt_mu_iso0p2",nbins,xmin,xmax);
  TH1F* num_SUSY_medPrompt_mu_iso4  = new TH1F("num_SUSY_medPrompt_mu_iso4","num_SUSY_medPrompt_mu_iso4",nbins,xmin,xmax);
  TH1F* num_SUSY_medPrompt_mu_iso6  = new TH1F("num_SUSY_medPrompt_mu_iso6","num_SUSY_medPrompt_mu_iso6",nbins,xmin,xmax);

  //Muon Isolation fakerate histograms
  TH1F* num_SUSY_medPrompt_mu_iso0p1_fr  = new TH1F("num_SUSY_medPrompt_mu_iso0p1_fr","num_SUSY_medPrompt_mu_iso0p1_fr",nbins,xmin,xmax);
  TH1F* num_SUSY_medPrompt_mu_iso0p2_fr  = new TH1F("num_SUSY_medPrompt_mu_iso0p2_fr","num_SUSY_medPrompt_mu_iso0p2_fr",nbins,xmin,xmax);
  TH1F* num_SUSY_medPrompt_mu_iso4_fr  = new TH1F("num_SUSY_medPrompt_mu_iso4_fr","num_SUSY_medPrompt_mu_iso4_fr",nbins,xmin,xmax);
  TH1F* num_SUSY_medPrompt_mu_iso6_fr  = new TH1F("num_SUSY_medPrompt_mu_iso6_fr","num_SUSY_medPrompt_mu_iso6_fr",nbins,xmin,xmax);

  //Denominator for electron and muon working points
  TH1F *den_tight_ele = new TH1F("den_tight_ele","den_tight_ele",nbins,xmin,xmax);
  TH1F *den_tight_fake_ele = new TH1F("den_tight_fake_ele","den_tight_fake_ele",nbins,xmin,xmax);
  TH1F *den_medPrompt_mu = new TH1F("den_medPrompt_mu","den_medPrompt_mu",nbins,xmin,xmax);
  TH1F *den_medPrompt_fake_mu = new TH1F("den_medPrompt_fake_mu","den_medPrompt_fake_mu",nbins,xmin,xmax);

  //Denominator isolation cut fake leptons
  TH1F *den_iso4_fake_mu = new TH1F("den_iso4_fake_mu","den_iso4_fake_mu",nbins,xmin,xmax);

  //Numerator isolation cut fake leptons
  TH1F *num_mp_iso4_fr_mu = new TH1F("num_mp_iso4_fr_mu","num_mp_iso4_fr_mu",nbins,xmin,xmax);
  TH1F *num_m_iso4_fr_mu = new TH1F("num_m_iso4_fr_mu","num_m_iso4_fr_mu",nbins,xmin,xmax);
  TH1F *num_t_iso4_fr_mu = new TH1F("num_t_iso4_fr_mu","num_t_iso4_fr_mu",nbins,xmin,xmax);
  TH1F *num_sid_iso4_fr_mu = new TH1F("num_sid_iso4_fr_mu","num_sid_iso4_fr_mu",nbins,xmin,xmax);
  TH1F *num_smva_iso4_fr_mu = new TH1F("num_smva_iso4_fr_mu","num_smva_iso4_fr_mu",nbins,xmin,xmax);

  string line;
  TChain *tch = new TChain("KUAnalysis");

  ifstream file (txtfile);
  if (file.is_open()){
    while (!file.eof()){
      getline(file,line);
      tch->Add(TString(line));
    }
  }

  KUAnalysis *ku = new KUAnalysis(tch);

  int nEntries = ku->fChain->GetEntries();
  if(test)
    nEntries = 100000;

  for(int e = 0; e < nEntries; e++){
    if (e % 10000 == 0) {
      fprintf(stdout, "\r  Processed events: %8d of %8d ", e, nEntries);
    }
    fflush(stdout);

    ku->fChain->GetEntry(e);

    weight = ku->weight;

    isGenZtoEle = false;
    isGenWtoEle = false;
    isGenZtoMu = false;
    isGenWtoMu = false;
    isGenTauToEle = false;
    isGenTauToMu = false;

    //Fill TlorentzVector with gen electrons that came from Z or W
    for(int g = 0; g < ku->genNele; g++){

	if(ku->genMomPDGID_ele->at(g) == 23)
	  isGenZtoEle = true;
	if(fabs(ku->genMomPDGID_ele->at(g)) == 24)
          isGenWtoEle = true;
	if(fabs(ku->genMomPDGID_ele->at(g)) == 15)
          isGenTauToEle = true;

	if(ku->genMomPDGID_ele->at(g) == 23 || fabs(ku->genMomPDGID_ele->at(g)) == 24){
	  gen_ele.SetPtEtaPhiM(ku->genPT_ele->at(g),ku->genEta_ele->at(g),ku->genPhi_ele->at(g),ku->genM_ele->at(g));
	  gen_ele_vec.push_back(gen_ele);
	}
    }

    //Fill TlorentzVector with gen muons that came from Z or W                                                                                            
    for(int g = 0; g < ku->genNmu; g++){

        if(ku->genMomPDGID_mu->at(g) == 23)
          isGenZtoMu = true;
        if(fabs(ku->genMomPDGID_mu->at(g)) == 24)
          isGenWtoMu = true;
	if(fabs(ku->genMomPDGID_mu->at(g)) == 15)
          isGenTauToMu = true;
 
	if(ku->genMomPDGID_mu->at(g) == 23 || fabs(ku->genMomPDGID_mu->at(g)) == 24){
	  gen_mu.SetPtEtaPhiM(ku->genPT_mu->at(g),ku->genEta_mu->at(g),ku->genPhi_mu->at(g),ku->genM_mu->at(g));
	  gen_mu_vec.push_back(gen_mu);
	}
    }

    isGenZ = isGenZtoEle || isGenZtoMu;
    isGenW = isGenWtoEle || isGenWtoMu;
    isGenTau = isGenTauToEle || isGenTauToMu;

    if(isGenZ || isGenW){

      //Electron gen matching
      if(isGenZtoEle || isGenWtoEle){
	for(int i = 0; i < ku->Nele; i++){
	  
	  ele_cand.SetPtEtaPhiM(ku->PT_ele->at(i),ku->Eta_ele->at(i),ku->Phi_ele->at(i),ku->M_ele->at(i));
	  for(int g = 0; g < gen_ele_vec.size(); g++){
	    dR_ele->Fill(ele_cand.DeltaR(gen_ele_vec[g]));
	    if(ele_cand.DeltaR(gen_ele_vec[g]) < dR_cut_ele){
	      passIdx_ele.push_back(i);
	      break;
	    }
	  }
	}
	gen_ele_vec.clear();
      }

      //Loop over gen matched electrons
      for(int i = 0; i < passIdx_ele.size(); i++){ 

	idx = passIdx_ele[i];

	pt = ku->Electron_pt->at(idx);
	miniIso = ku->MiniIso_ele->at(idx);

	miniIsoVsPt_sig->Fill(pt,miniIso,weight);
	miniIsoPtVsPt_sig->Fill(pt,miniIso*pt,weight);
	//cout <<  ku->weight << endl; 

	if(ku->ID_ele->at(idx) >= 0 && pt > 5){
	  den_all_ele->Fill(pt,weight);
	  //SUSY wp
	  if(ku->FO_baseline_ele->at(idx) && ku->baseline_loose_ele->at(idx))
	    num_vl_ele->Fill(pt,weight);
	  if(ku->ID_ele->at(idx) >= 2)
	    num_l_ele->Fill(pt,weight);
	  if(ku->ID_ele->at(idx) >= 3)
	    num_m_ele->Fill(pt,weight);
	  if(ku->ID_ele->at(idx) >= 4)// && ku->baseline_tight_ele->at(idx))
	    num_t_ele->Fill(pt,weight);

	  //SUSY wp with baseline criteria
	  if(ku->FO_baseline_ele->at(idx) && ku->baseline_loose_ele->at(idx)){
	    if(ku->ID_ele->at(idx) >= 2)
	      num_vl_bl_ele->Fill(pt,weight);
	    if(ku->ID_ele->at(idx) >= 3)
	      num_l_bl_ele->Fill(pt,weight);
	    if(ku->ID_ele->at(idx) >= 4)// && ku->baseline_tight_ele->at(idx))
	      num_t_bl_ele->Fill(pt,weight);
	  }

	  //EGamma wp Fall17V1
	  if(ku->Electron_mvaFall17V1noIso_WP80->at(idx))
            num_WP80V1->Fill(pt,weight);
          if(ku->Electron_mvaFall17V1noIso_WP90->at(idx))
            num_WP90V1->Fill(pt,weight);
          if(ku->Electron_mvaFall17V1noIso_WPL->at(idx))
            num_WPLV1->Fill(pt,weight);
	  //EGamma wp Fall17V2
          if(ku->Electron_mvaFall17V2noIso_WP80->at(idx))
            num_WP80V2->Fill(pt,weight);
          if(ku->Electron_mvaFall17V2noIso_WP90->at(idx))
            num_WP90V2->Fill(pt,weight);
          if(ku->Electron_mvaFall17V2noIso_WPL->at(idx))
            num_WPLV2->Fill(pt,weight);

	  //electron candidate requirements
	  if(ku->Electron_lostHits->at(idx) == 0)
            num_lostHits->Fill(pt,weight);
	  if(ku->Electron_ip3d->at(idx) < ip3D)
	    num_ip3D->Fill(pt,weight);
	  if(ku->Electron_sip3d->at(idx) < sip3D)
	    num_sip3D->Fill(pt,weight);
	  if(ku->Electron_sip3d->at(idx) < 4.0)
            num_sip3D_4->Fill(pt,weight);
	  if(ku->Electron_sip3d->at(idx) < 8.0)
            num_sip3D_8->Fill(pt,weight);
	  if(ku->Electron_convVeto->at(idx))
	    num_convVeto->Fill(pt,weight);
	  if(fabs(ku->Electron_dxy->at(idx)) < dxy)
	    num_dxy->Fill(pt,weight);
	  if(fabs(ku->Electron_dz->at(idx)) < dz)
	    num_dz->Fill(pt,weight);
	}

	//Tight electron denominator
	if(ku->ID_ele->at(idx) >= 4 && pt > 5){
	  den_tight_ele->Fill(pt,weight);
	  if(miniIso < 0.1)
	    num_SUSY_tight_ele_iso0p1->Fill(pt,weight);
	  if(miniIso < 0.2)
	    num_SUSY_tight_ele_iso0p2->Fill(pt,weight);
	  if(miniIso*pt < 4.)
	    num_SUSY_tight_ele_iso4->Fill(pt,weight);
	  if(miniIso*pt < 6.)
	    num_SUSY_tight_ele_iso6->Fill(pt,weight);
	}
      }
      passIdx_ele.clear();      

      //Muon gen matching                                                                                                                                 
      if(isGenZtoMu ||isGenWtoMu){
	for(int i = 0; i < ku->Nmu; i++){
	
	  mu_cand.SetPtEtaPhiM(ku->PT_mu->at(i),ku->Eta_mu->at(i),ku->Phi_mu->at(i),ku->M_mu->at(i));
	  for(int g = 0; g < gen_mu_vec.size(); g++){
	    dR_mu->Fill(mu_cand.DeltaR(gen_mu_vec[g]));
	    if(mu_cand.DeltaR(gen_mu_vec[g]) < dR_cut_mu){
	      passIdx_mu.push_back(i);
	      break;
	    }
	  }
	}
	gen_mu_vec.clear();
      }

      //Loop over gen matched muons
      for(int i = 0; i < passIdx_mu.size(); i++){  
	
	idx = passIdx_mu[i];

	pt = ku->Muon_pt->at(idx);
        miniIso = ku->MiniIso_mu->at(idx);

	miniIsoVsPt_mu_sig->Fill(pt,miniIso,weight);
        miniIsoPtVsPt_mu_sig->Fill(pt,miniIso*pt,weight);

	if(ku->ID_mu->at(idx) >= 0 && pt > 5){
          den_all_mu->Fill(pt,weight);
	  //Muon working points 
	  if(ku->baseline_loose_mu->at(idx))
	    num_bl_mu->Fill(pt,weight);
	  if(ku->Muon_mediumId->at(idx))
	    num_m_mu->Fill(pt,weight);
	  if(ku->Muon_mediumPromptId->at(idx))
            num_mp_mu->Fill(pt,weight);
	  if(ku->Muon_tightId->at(idx))// && ku->baseline_tight_mu->at(idx))
            num_t_mu->Fill(pt,weight);
	  if(ku->Muon_softId->at(idx))
            num_sid_mu->Fill(pt,weight);
	  if(ku->Muon_softMvaId->at(idx))
            num_smva_mu->Fill(pt,weight);
	  
	  //Muon selection criteria
          if(ku->Muon_ip3d->at(idx) < ip3D)
            num_ip3D_mu->Fill(pt,weight);
          if(ku->Muon_sip3d->at(idx) < sip3D)
            num_sip3D_mu->Fill(pt,weight);
	  if(ku->Muon_sip3d->at(idx) < 4.0)
            num_sip3D_mu_4->Fill(pt,weight);
	  if(ku->Muon_sip3d->at(idx) < 8.0)
            num_sip3D_mu_8->Fill(pt,weight);
          if(fabs(ku->Muon_dxy->at(idx)) < dxy)
            num_dxy_mu->Fill(pt,weight);
          if(fabs(ku->Muon_dz->at(idx)) < dz)
            num_dz_mu->Fill(pt,weight);
	}
	if(ku->Muon_mediumPromptId->at(idx) && pt > 5){
          den_medPrompt_mu->Fill(pt,weight);
          if(miniIso < 0.1)
            num_SUSY_medPrompt_mu_iso0p1->Fill(pt,weight);
          if(miniIso < 0.2)
            num_SUSY_medPrompt_mu_iso0p2->Fill(pt,weight);
          if(miniIso*pt < 4.)
            num_SUSY_medPrompt_mu_iso4->Fill(pt,weight);
          if(miniIso*pt < 6.)
            num_SUSY_medPrompt_mu_iso6->Fill(pt,weight);
        }
      }
      passIdx_mu.clear();
      
    }
    if(!isGenZ && !isGenW && !isGenTau){
      //cout << "electron fake ID's" << endl;
      //Electron fakes
      //if(ku->genNele == 0){//For DYJets only
	for(int i = 0; i < ku->Nele; i++){

	  pt = ku->Electron_pt->at(i);
	  miniIso = ku->MiniIso_ele->at(i);
	  
	  miniIsoVsPt_fake->Fill(pt,miniIso,weight);
	  miniIsoPtVsPt_fake->Fill(pt,miniIso*pt,weight);
	  
	  if(ku->ID_ele->at(i) >= 0 && pt > 5){
	    //for(int j = 0; j < ku->genNele; j++)
	      //cout << "electron from: " << ku->genMomPDGID_ele->at(j) << endl;
	    den_fake_ele->Fill(pt,weight);
	    //Selection criteria fakerates
	    if(ku->Electron_lostHits->at(i) == 0)
	      num_lostHits_f->Fill(pt,weight);
	    if(ku->Electron_ip3d->at(i) < ip3D)
	      num_ip3D_f->Fill(pt,weight);
	    if(ku->Electron_sip3d->at(i) < sip3D)
	      num_sip3D_f->Fill(pt,weight);
	    if(ku->Electron_sip3d->at(i) < 4.0)
	      num_sip3D_f_4->Fill(pt,weight);
	    if(ku->Electron_sip3d->at(i) < 8.0)
	      num_sip3D_f_8->Fill(pt,weight);
	    if(ku->Electron_convVeto->at(i))
	      num_convVeto_f->Fill(pt,weight);
	    if(fabs(ku->Electron_dxy->at(i)) < dxy)
	      num_dxy_f->Fill(pt,weight);
	    if(fabs(ku->Electron_dz->at(i)) < dz)
	      num_dz_f->Fill(pt,weight);
	    
	    //SUSY wp fakerates
	    if(ku->FO_baseline_ele->at(i) && ku->baseline_loose_ele->at(i))
	      num_bl_fr_ele->Fill(pt,weight);
	    if(ku->ID_ele->at(i) >= 2)
	      num_vl_fr_ele->Fill(pt,weight);
	    if(ku->ID_ele->at(i) >= 3)
	      num_l_fr_ele->Fill(pt,weight);
	    if(ku->ID_ele->at(i) >= 4)// && ku->baseline_tight_ele->at(i))
	      num_t_fr_ele->Fill(pt,weight);
	    
	    //EGamma wp Fall17V1 fakerates
	    if(ku->Electron_mvaFall17V1noIso_WP80->at(i))
	      num_WP80V1_fr->Fill(pt,weight);
	    if(ku->Electron_mvaFall17V1noIso_WP90->at(i))
	      num_WP90V1_fr->Fill(pt,weight);
	    if(ku->Electron_mvaFall17V1noIso_WPL->at(i))
	      num_WPLV1_fr->Fill(pt,weight);
	    
	    //EGamma wp Fall17V1 fakerates
	    if(ku->Electron_mvaFall17V2noIso_WP80->at(i))
	      num_WP80V2_fr->Fill(pt,weight);
	    if(ku->Electron_mvaFall17V2noIso_WP90->at(i))
	      num_WP90V2_fr->Fill(pt,weight);
	    if(ku->Electron_mvaFall17V2noIso_WPL->at(i))
	      num_WPLV2_fr->Fill(pt,weight);
	  }
	  if(ku->ID_ele->at(i) >= 4 && pt > 5){
	    den_tight_fake_ele->Fill(pt,weight);
	    if(miniIso < 0.1)
	      num_SUSY_tight_ele_iso0p1_fr->Fill(pt,weight);
	    if(miniIso < 0.2)
	      num_SUSY_tight_ele_iso0p2_fr->Fill(pt,weight);
	    if(miniIso*pt < 4.)
	      num_SUSY_tight_ele_iso4_fr->Fill(pt,weight);
	    if(miniIso*pt < 6.)
	      num_SUSY_tight_ele_iso6_fr->Fill(pt,weight);
	  }
	  //}
      }//end electron fakes
      
      //Muon fakes
	if(ku->genNmu == 0){//For TTJets and DYJets
	//cout << "muon fake ID's" << endl;
	for(int i = 0; i < ku->Nmu; i++){

	  pt = ku->Muon_pt->at(i);
	  miniIso = ku->MiniIso_mu->at(i);

	  miniIsoVsPt_mu_fake->Fill(pt,miniIso,weight);
	  miniIsoPtVsPt_mu_fake->Fill(pt,miniIso*pt,weight);
	  
	  if(ku->ID_mu->at(i) >= 0 && pt > 5){
	    //for(int j = 0; j < ku->genNmu; j++)
	      //cout << "muon from: " << ku->genMomPDGID_mu->at(j) << endl;
	    den_fake_mu->Fill(pt,weight);
	    if(ku->baseline_loose_mu->at(i))
	      num_bl_fr_mu->Fill(pt,weight);
	    if(ku->Muon_mediumId->at(i))
	      num_m_fr_mu->Fill(pt,weight);
	    if(ku->Muon_mediumPromptId->at(i))
	      num_mp_fr_mu->Fill(pt,weight);
	    if(ku->Muon_tightId->at(i))// && ku->baseline_tight_mu->at(i))
	      num_t_fr_mu->Fill(pt,weight);
	    if(ku->Muon_softId->at(i))
	      num_sid_fr_mu->Fill(pt,weight);
	    if(ku->Muon_softMvaId->at(i))
	      num_smva_fr_mu->Fill(pt,weight);

	    //Muon selection criteria
	    if(ku->Muon_ip3d->at(i) < ip3D)
	      num_ip3D_fr_mu->Fill(pt,weight);
	    if(ku->Muon_sip3d->at(i) < sip3D)
	      num_sip3D_fr_mu->Fill(pt,weight);
	    if(ku->Muon_sip3d->at(i) < 4.0)
	      num_sip3D_fr_mu_4->Fill(pt,weight);
	    if(ku->Muon_sip3d->at(i) < 8.0)
	      num_sip3D_fr_mu_8->Fill(pt,weight);
	    if(fabs(ku->Muon_dxy->at(i)) < dxy)
	      num_dxy_fr_mu->Fill(pt,weight);
	    if(fabs(ku->Muon_dz->at(i)) < dz)
	      num_dz_fr_mu->Fill(pt,weight);
	  }
	  if(ku->Muon_mediumPromptId->at(i) && pt > 5){
	    den_medPrompt_fake_mu->Fill(pt,weight);
	    if(miniIso < 0.1)
	      num_SUSY_medPrompt_mu_iso0p1_fr->Fill(pt,weight);
	    if(miniIso < 0.2)
	      num_SUSY_medPrompt_mu_iso0p2_fr->Fill(pt,weight);
	    if(miniIso*pt < 4.)
	      num_SUSY_medPrompt_mu_iso4_fr->Fill(pt,weight);
	    if(miniIso*pt < 6.)
	      num_SUSY_medPrompt_mu_iso6_fr->Fill(pt,weight);
	  }

	  if(miniIso*pt < 4. && pt > 5){
	    den_iso4_fake_mu->Fill(pt,weight);
	    if(ku->Muon_mediumId->at(i))
	      num_m_iso4_fr_mu->Fill(pt,weight);
	    if(ku->Muon_mediumPromptId->at(i))
	      num_mp_iso4_fr_mu->Fill(pt,weight);
	    if(ku->Muon_tightId->at(i))
	      num_t_iso4_fr_mu->Fill(pt,weight);
	    if(ku->Muon_softId->at(i))
	      num_sid_iso4_fr_mu->Fill(pt,weight);
	    if(ku->Muon_softMvaId->at(i))
	      num_smva_iso4_fr_mu->Fill(pt,weight);
	  }
	}
      }//end muon fakes

    }
  }
  cout << endl;

  //Working point efficiencies (SUSY group)
  vector<TEfficiency*> eff_wp_ele;
  vector<TString> label_wp_ele;

  TEfficiency *eff_vl_ele = new TEfficiency(*num_vl_ele,*den_all_ele);
  TEfficiency *eff_l_ele = new TEfficiency(*num_l_ele,*den_all_ele);
  TEfficiency *eff_m_ele = new TEfficiency(*num_m_ele,*den_all_ele);
  TEfficiency *eff_t_ele = new TEfficiency(*num_t_ele,*den_all_ele);

  //eff_wp_ele.push_back(eff_vl_ele); label_wp_ele.push_back("baseline");
  eff_wp_ele.push_back(eff_l_ele); label_wp_ele.push_back("very loose");
  eff_wp_ele.push_back(eff_m_ele); label_wp_ele.push_back("loose");
  eff_wp_ele.push_back(eff_t_ele); label_wp_ele.push_back("tight");

  PlotEff(eff_wp_ele,label_wp_ele,"eff_wp_ele", ele_label, ylabel_eff, samplename, true);
  PlotEff(eff_wp_ele,label_wp_ele,"eff_wp_ele_zoom", ele_label, ylabel_eff, samplename, true, 1.01, 0.7);

  //Working point fakerates (SUSY group)
  vector<TEfficiency*> fr_wp_ele;

  TEfficiency *fr_bl_ele = new TEfficiency(*num_bl_fr_ele,*den_fake_ele);
  TEfficiency *fr_vl_ele = new TEfficiency(*num_vl_fr_ele,*den_fake_ele);
  TEfficiency *fr_l_ele = new TEfficiency(*num_l_fr_ele,*den_fake_ele);
  TEfficiency *fr_t_ele = new TEfficiency(*num_t_fr_ele,*den_fake_ele);

  //fr_wp_ele.push_back(fr_bl_ele); 
  fr_wp_ele.push_back(fr_vl_ele); 
  fr_wp_ele.push_back(fr_l_ele); 
  fr_wp_ele.push_back(fr_t_ele); 

  PlotEff(fr_wp_ele,label_wp_ele,"fr_wp_ele",ele_label, ylabel_fr, samplename, false);
  PlotEff(fr_wp_ele,label_wp_ele,"fr_wp_ele_zoom",ele_label, ylabel_fr, samplename, false, 0.4);

  //Working point efficiencies with baseline selection (SUSY group)
  vector<TEfficiency*> eff_wp_bl_ele;
  vector<TString> label_wp_bl_ele;

  TEfficiency *eff_vl_bl_ele = new TEfficiency(*num_vl_bl_ele,*den_all_ele);
  TEfficiency *eff_l_bl_ele = new TEfficiency(*num_l_bl_ele,*den_all_ele);
  TEfficiency *eff_t_bl_ele = new TEfficiency(*num_t_bl_ele,*den_all_ele);

  eff_wp_bl_ele.push_back(eff_vl_bl_ele); label_wp_bl_ele.push_back("very loose");
  eff_wp_bl_ele.push_back(eff_l_bl_ele); label_wp_bl_ele.push_back("loose");
  eff_wp_bl_ele.push_back(eff_t_bl_ele); label_wp_bl_ele.push_back("tight");

  //PlotEff(eff_wp_bl_ele,label_wp_bl_ele,"eff_wp_bl_ele", ele_label,ylabel_eff, samplename, true);

  //Working point efficiencies (Egamma Fall17V1)
  vector<TEfficiency*> eff_wp_v1;
  vector<TString> label_wp_egamma;

  TEfficiency *eff_wp80_v1 = new TEfficiency(*num_WP80V1,*den_all_ele);
  TEfficiency *eff_wp90_v1 = new TEfficiency(*num_WP90V1,*den_all_ele);
  TEfficiency *eff_wpl_v1 = new TEfficiency(*num_WPLV1,*den_all_ele);

  eff_wp_v1.push_back(eff_wpl_v1); label_wp_egamma.push_back("WPL");  
  eff_wp_v1.push_back(eff_wp90_v1); label_wp_egamma.push_back("WP90");
  eff_wp_v1.push_back(eff_wp80_v1); label_wp_egamma.push_back("WP80");

  PlotEff(eff_wp_v1,label_wp_egamma,"eff_wp_egV1", ele_label,ylabel_eff, samplename, true);

  //Working point fakerate (Egamma Fall17V1)
  vector<TEfficiency*> fr_wp_v1;

  TEfficiency *fr_wp80_v1 = new TEfficiency(*num_WP80V1_fr,*den_fake_ele);
  TEfficiency *fr_wp90_v1 = new TEfficiency(*num_WP90V1_fr,*den_fake_ele);
  TEfficiency *fr_wpl_v1 = new TEfficiency(*num_WPLV1_fr,*den_fake_ele);

  fr_wp_v1.push_back(fr_wpl_v1); 
  fr_wp_v1.push_back(fr_wp90_v1); 
  fr_wp_v1.push_back(fr_wp80_v1); 

  PlotEff(fr_wp_v1,label_wp_egamma,"fr_wp_egV1", ele_label,ylabel_fr, samplename, false);

  //Working point efficiencies (Egamma Fall17V2)
  vector<TEfficiency*> eff_wp_v2;

  TEfficiency *eff_wp80_v2 = new TEfficiency(*num_WP80V2,*den_all_ele);
  TEfficiency *eff_wp90_v2 = new TEfficiency(*num_WP90V2,*den_all_ele);
  TEfficiency *eff_wpl_v2 = new TEfficiency(*num_WPLV2,*den_all_ele);

  eff_wp_v2.push_back(eff_wpl_v2);
  eff_wp_v2.push_back(eff_wp90_v2); 
  eff_wp_v2.push_back(eff_wp80_v2);

  //PlotEff(eff_wp_v2,label_wp_egamma,"eff_wp_egV2", ele_label,ylabel_eff, samplename, true);

  //Working point fakerate (Egamma Fall17V2)
  vector<TEfficiency*> fr_wp_v2;

  TEfficiency *fr_wp80_v2 = new TEfficiency(*num_WP80V2_fr,*den_fake_ele);
  TEfficiency *fr_wp90_v2 = new TEfficiency(*num_WP90V2_fr,*den_fake_ele);
  TEfficiency *fr_wpl_v2 = new TEfficiency(*num_WPLV2_fr,*den_fake_ele);

  fr_wp_v2.push_back(fr_wpl_v2);
  fr_wp_v2.push_back(fr_wp90_v2);
  fr_wp_v2.push_back(fr_wp80_v2);

  //PlotEff(fr_wp_v2,label_wp_egamma,"fr_wp_egV2", ele_label, ylabel_fr, samplename, false);

  //Working point muon efficiencies
  vector<TEfficiency*> eff_wp_mu;
  vector<TString> label_wp_mu;

  TEfficiency *eff_bl_mu = new TEfficiency(*num_bl_mu,*den_all_mu);
  TEfficiency *eff_m_mu = new TEfficiency(*num_m_mu,*den_all_mu);
  TEfficiency *eff_mp_mu = new TEfficiency(*num_mp_mu,*den_all_mu);
  TEfficiency *eff_t_mu = new TEfficiency(*num_t_mu,*den_all_mu);
  TEfficiency *eff_sid_mu = new TEfficiency(*num_sid_mu,*den_all_mu);
  TEfficiency *eff_smva_mu = new TEfficiency(*num_smva_mu,*den_all_mu);

  //eff_wp_mu.push_back(eff_bl_mu); label_wp_mu.push_back("baseline");
  eff_wp_mu.push_back(eff_m_mu); label_wp_mu.push_back("medium");
  eff_wp_mu.push_back(eff_mp_mu); label_wp_mu.push_back("medium prompt");
  eff_wp_mu.push_back(eff_t_mu); label_wp_mu.push_back("tight");
  //eff_wp_mu.push_back(eff_sid_mu); label_wp_mu.push_back("soft ID");
  //eff_wp_mu.push_back(eff_smva_mu); label_wp_mu.push_back("soft MVA ID");

  PlotEff(eff_wp_mu, label_wp_mu, "eff_wp_mu", mu_label, ylabel_eff, samplename, true);
  PlotEff(eff_wp_mu, label_wp_mu, "eff_wp_mu_zoom", mu_label, ylabel_eff, samplename, true, 1.005, 0.9);

  //Working point muon fakerate
  vector<TEfficiency*> fr_wp_mu;

  TEfficiency *fr_bl_mu = new TEfficiency(*num_bl_fr_mu,*den_fake_mu);
  TEfficiency *fr_m_mu = new TEfficiency(*num_m_fr_mu,*den_fake_mu);
  TEfficiency *fr_mp_mu = new TEfficiency(*num_mp_fr_mu,*den_fake_mu);
  TEfficiency *fr_t_mu = new TEfficiency(*num_t_fr_mu,*den_fake_mu);
  TEfficiency *fr_sid_mu = new TEfficiency(*num_sid_fr_mu,*den_fake_mu);
  TEfficiency *fr_smva_mu = new TEfficiency(*num_smva_fr_mu,*den_fake_mu);

  //fr_wp_mu.push_back(fr_bl_mu); 
  fr_wp_mu.push_back(fr_m_mu); 
  fr_wp_mu.push_back(fr_mp_mu);
  fr_wp_mu.push_back(fr_t_mu); 
  //fr_wp_mu.push_back(fr_sid_mu);
  //fr_wp_mu.push_back(fr_smva_mu);

  PlotEff(fr_wp_mu, label_wp_mu, "fr_wp_mu", mu_label, ylabel_fr, samplename, false);
  PlotEff(fr_wp_mu, label_wp_mu, "fr_wp_mu_zoom", mu_label, ylabel_fr, samplename, false, 0.4);

  //Fakerate isolated muons (pt*miniIso < 4.0)
  vector<TEfficiency*> fr_iso4_wp_mu;
  vector<TString> label_iso4_wp_fr_mu;

  TEfficiency *fr_iso4_m_mu = new TEfficiency(*num_m_iso4_fr_mu,*den_iso4_fake_mu);
  TEfficiency *fr_iso4_mp_mu = new TEfficiency(*num_mp_iso4_fr_mu,*den_iso4_fake_mu);
  TEfficiency *fr_iso4_t_mu = new TEfficiency(*num_t_iso4_fr_mu,*den_iso4_fake_mu);
  TEfficiency *fr_iso4_sid_mu = new TEfficiency(*num_sid_iso4_fr_mu,*den_iso4_fake_mu);
  TEfficiency *fr_iso4_smva_mu = new TEfficiency(*num_smva_iso4_fr_mu,*den_iso4_fake_mu);

  fr_iso4_wp_mu.push_back(fr_iso4_m_mu);label_iso4_wp_fr_mu.push_back("medium");
  fr_iso4_wp_mu.push_back(fr_iso4_mp_mu);label_iso4_wp_fr_mu.push_back("medium prompt");
  fr_iso4_wp_mu.push_back(fr_iso4_t_mu);label_iso4_wp_fr_mu.push_back("tight");
  //fr_iso4_wp_mu.push_back(fr_iso4_sid_mu);label_iso4_wp_fr_mu.push_back("soft ID");
  //fr_iso4_wp_mu.push_back(fr_iso4_smva_mu);label_iso4_wp_fr_mu.push_back("soft MVA ID");

  PlotEff(fr_iso4_wp_mu, label_iso4_wp_fr_mu, "fr_iso4_wp_fr_mu", mu_label, ylabel_fr, samplename, true);

  //electron candidate selection efficiencies and fakerates
  vector<TEfficiency*> eff_sel;
  vector<TEfficiency*> fake_sel;
  vector<TString> label_sel;

  TEfficiency *eff_lostHits = new TEfficiency(*num_lostHits,*den_all_ele);
  TEfficiency *eff_ip3D = new TEfficiency(*num_ip3D,*den_all_ele);
  TEfficiency *eff_sip3D = new TEfficiency(*num_sip3D,*den_all_ele);
  TEfficiency *eff_sip3D_4 = new TEfficiency(*num_sip3D_4,*den_all_ele);
  TEfficiency *eff_sip3D_8 = new TEfficiency(*num_sip3D_8,*den_all_ele);
  TEfficiency *eff_convVeto = new TEfficiency(*num_convVeto,*den_all_ele);
  TEfficiency *eff_dxy = new TEfficiency(*num_dxy,*den_all_ele);
  TEfficiency *eff_dz = new TEfficiency(*num_dz,*den_all_ele);

  //eff_sel.push_back(eff_ip3D); label_sel.push_back(Form("ip3D < %.4f",ip3D));
  eff_sel.push_back(eff_sip3D); label_sel.push_back(Form("sip3D < %.1f",sip3D));
  //eff_sel.push_back(eff_sip3D_4); label_sel.push_back("sip3D < 4.0");
  //eff_sel.push_back(eff_sip3D_8); label_sel.push_back("sip3D < 8.0");
  eff_sel.push_back(eff_dxy); label_sel.push_back(Form("dxy < %.2f",dxy));
  eff_sel.push_back(eff_dz); label_sel.push_back(Form("dz < %.2f",dz));
  eff_sel.push_back(eff_lostHits); label_sel.push_back(Form("pass lostHits"));
  eff_sel.push_back(eff_convVeto); label_sel.push_back("pass convVeto");

  TEfficiency *fake_lostHits = new TEfficiency(*num_lostHits_f,*den_fake_ele);
  TEfficiency *fake_ip3D = new TEfficiency(*num_ip3D_f,*den_fake_ele);
  TEfficiency *fake_sip3D = new TEfficiency(*num_sip3D_f,*den_fake_ele);
  TEfficiency *fake_sip3D_4 = new TEfficiency(*num_sip3D_f_4,*den_fake_ele);
  TEfficiency *fake_sip3D_8 = new TEfficiency(*num_sip3D_f_8,*den_fake_ele);
  TEfficiency *fake_convVeto = new TEfficiency(*num_convVeto_f,*den_fake_ele);
  TEfficiency *fake_dxy = new TEfficiency(*num_dxy_f,*den_fake_ele);
  TEfficiency *fake_dz = new TEfficiency(*num_dz_f,*den_fake_ele);

  //fake_sel.push_back(fake_ip3D);
  fake_sel.push_back(fake_sip3D);
  //fake_sel.push_back(fake_sip3D_4);
  //fake_sel.push_back(fake_sip3D_8);
  fake_sel.push_back(fake_dxy);
  fake_sel.push_back(fake_dz);
  fake_sel.push_back(fake_lostHits);
  fake_sel.push_back(fake_convVeto);

  PlotEff(eff_sel,label_sel,"eff_sel", ele_label,ylabel_eff, samplename, true);
  PlotEff(fake_sel,label_sel,"fake_sel", ele_label,ylabel_fr, samplename, true);

  //muon candidate selection efficiencies and fakerates
  vector<TEfficiency*> eff_sel_mu;
  vector<TEfficiency*> fake_sel_mu;
  vector<TString> label_sel_mu;

  TEfficiency *eff_ip3D_mu  = new TEfficiency(*num_ip3D_mu,*den_all_mu);
  TEfficiency *eff_sip3D_mu = new TEfficiency(*num_sip3D_mu,*den_all_mu);
  TEfficiency *eff_sip3D_mu_4 = new TEfficiency(*num_sip3D_mu_4,*den_all_mu);
  TEfficiency *eff_sip3D_mu_8 = new TEfficiency(*num_sip3D_mu_8,*den_all_mu);
  TEfficiency *eff_dxy_mu   = new TEfficiency(*num_dxy_mu,*den_all_mu);
  TEfficiency *eff_dz_mu    = new TEfficiency(*num_dz_mu,*den_all_mu);

  //eff_sel_mu.push_back(eff_ip3D_mu); label_sel_mu.push_back(Form("ip3D < %.4f",ip3D));
  eff_sel_mu.push_back(eff_sip3D_mu); label_sel_mu.push_back(Form("sip3D < %.1f",sip3D));
  //eff_sel_mu.push_back(eff_sip3D_mu_4); label_sel_mu.push_back("sip3D < 4.0");
  //eff_sel_mu.push_back(eff_sip3D_mu_8); label_sel_mu.push_back("sip3D < 8.0");
  eff_sel_mu.push_back(eff_dxy_mu); label_sel_mu.push_back(Form("dxy < %.2f",dxy));
  eff_sel_mu.push_back(eff_dz_mu); label_sel_mu.push_back(Form("dz < %.2f",dz));

  TEfficiency *fake_ip3D_mu  = new TEfficiency(*num_ip3D_fr_mu,*den_fake_mu);
  TEfficiency *fake_sip3D_mu = new TEfficiency(*num_sip3D_fr_mu,*den_fake_mu);
  TEfficiency *fake_sip3D_mu_4 = new TEfficiency(*num_sip3D_fr_mu_4,*den_fake_mu);
  TEfficiency *fake_sip3D_mu_8 = new TEfficiency(*num_sip3D_fr_mu_8,*den_fake_mu);
  TEfficiency *fake_dxy_mu   = new TEfficiency(*num_dxy_fr_mu,*den_fake_mu);
  TEfficiency *fake_dz_mu    = new TEfficiency(*num_dz_fr_mu,*den_fake_mu);

  //fake_sel_mu.push_back(fake_ip3D_mu);
  fake_sel_mu.push_back(fake_sip3D_mu);
  //fake_sel_mu.push_back(fake_sip3D_mu_4);
  //fake_sel_mu.push_back(fake_sip3D_mu_8);
  fake_sel_mu.push_back(fake_dxy_mu);
  fake_sel_mu.push_back(fake_dz_mu);

  PlotEff(eff_sel_mu,label_sel_mu,"eff_sel_mu", mu_label,ylabel_eff, samplename, true);
  PlotEff(fake_sel_mu,label_sel_mu,"fake_sel_mu", mu_label,ylabel_fr, samplename, true);

  //Electron WP comparison (SUSY and Fall17V2)
  vector<TEfficiency*> eff_wp_comb_v2_ele;
  vector<TString> label_wp_comb_v2_ele;

  eff_wp_comb_v2_ele.push_back(eff_l_ele); label_wp_comb_v2_ele.push_back("very loose");
  eff_wp_comb_v2_ele.push_back(eff_m_ele); label_wp_comb_v2_ele.push_back("loose");
  eff_wp_comb_v2_ele.push_back(eff_t_ele); label_wp_comb_v2_ele.push_back("tight");
  eff_wp_comb_v2_ele.push_back(eff_wpl_v2); label_wp_comb_v2_ele.push_back("WPL");
  eff_wp_comb_v2_ele.push_back(eff_wp80_v2); label_wp_comb_v2_ele.push_back("WP80");
  eff_wp_comb_v2_ele.push_back(eff_wp90_v2); label_wp_comb_v2_ele.push_back("WP90");

  PlotEff(eff_wp_comb_v2_ele,label_wp_comb_v2_ele,"eff_wp_comb_v2_ele", ele_label, ylabel_eff, samplename, true);
  PlotEff(eff_wp_comb_v2_ele,label_wp_comb_v2_ele,"eff_wp_comb_v2_ele_zoom", ele_label, ylabel_eff, samplename, true, 1.01, 0.5);

  //Electron WP comparison fakerate (SUSY and Fall17V2)
  vector<TEfficiency*> fr_wp_comb_v2_ele;

  fr_wp_comb_v2_ele.push_back(fr_vl_ele); 
  fr_wp_comb_v2_ele.push_back(fr_l_ele); 
  fr_wp_comb_v2_ele.push_back(fr_t_ele); 
  fr_wp_comb_v2_ele.push_back(fr_wpl_v2); 
  fr_wp_comb_v2_ele.push_back(fr_wp80_v2);
  fr_wp_comb_v2_ele.push_back(fr_wp90_v2);

  PlotEff(fr_wp_comb_v2_ele,label_wp_comb_v2_ele,"fr_wp_comb_v2_ele", ele_label, ylabel_fr, samplename, false);
  PlotEff(fr_wp_comb_v2_ele,label_wp_comb_v2_ele,"fr_wp_comb_v2_ele_zoom", ele_label, ylabel_fr, samplename, false,0.4);

  //Electron tight SUSY WP isolation efficiency
  vector<TEfficiency*> eff_tight_ele_iso;
  vector<TString> label_tight_ele_iso;

  TEfficiency *eff_tight_ele_iso0p1  = new TEfficiency(*num_SUSY_tight_ele_iso0p1,*den_tight_ele);
  TEfficiency *eff_tight_ele_iso0p2  = new TEfficiency(*num_SUSY_tight_ele_iso0p2,*den_tight_ele);
  TEfficiency *eff_tight_ele_iso4  = new TEfficiency(*num_SUSY_tight_ele_iso4,*den_tight_ele);
  TEfficiency *eff_tight_ele_iso6  = new TEfficiency(*num_SUSY_tight_ele_iso6,*den_tight_ele);

  //eff_tight_ele_iso.push_back(eff_t_ele); label_tight_ele_iso.push_back("tight (no iso)");
  eff_tight_ele_iso.push_back(eff_tight_ele_iso0p1); label_tight_ele_iso.push_back("miniIso < 0.1");
  eff_tight_ele_iso.push_back(eff_tight_ele_iso0p2); label_tight_ele_iso.push_back("miniIso < 0.2");
  eff_tight_ele_iso.push_back(eff_tight_ele_iso4); label_tight_ele_iso.push_back("miniIso #upoint p_{T} < 4");
  eff_tight_ele_iso.push_back(eff_tight_ele_iso6); label_tight_ele_iso.push_back("miniIso #upoint p_{T} < 6");

  PlotEff(eff_tight_ele_iso,label_tight_ele_iso,"eff_tight_ele_iso", ele_label, ylabel_eff, samplename, true);

  //Electron tight SUSY WP isolation fakerate
  vector<TEfficiency*> fr_tight_ele_iso;

  TEfficiency *fr_tight_ele_iso0p1  = new TEfficiency(*num_SUSY_tight_ele_iso0p1_fr,*den_tight_fake_ele);
  TEfficiency *fr_tight_ele_iso0p2  = new TEfficiency(*num_SUSY_tight_ele_iso0p2_fr,*den_tight_fake_ele);
  TEfficiency *fr_tight_ele_iso4  = new TEfficiency(*num_SUSY_tight_ele_iso4_fr,*den_tight_fake_ele);
  TEfficiency *fr_tight_ele_iso6  = new TEfficiency(*num_SUSY_tight_ele_iso6_fr,*den_tight_fake_ele);

  //fr_tight_ele_iso.push_back(fr_t_ele);
  fr_tight_ele_iso.push_back(fr_tight_ele_iso0p1); 
  fr_tight_ele_iso.push_back(fr_tight_ele_iso0p2); 
  fr_tight_ele_iso.push_back(fr_tight_ele_iso4); 
  fr_tight_ele_iso.push_back(fr_tight_ele_iso6); 

  PlotEff(fr_tight_ele_iso,label_tight_ele_iso,"fr_tight_ele_iso", ele_label, ylabel_fr, samplename, false);

  //Muon medium prompt WP isolation efficiency
  vector<TEfficiency*> eff_medPrompt_mu_iso;
  vector<TString> label_medPrompt_mu_iso;

  TEfficiency *eff_medPrompt_mu_iso0p1  = new TEfficiency(*num_SUSY_medPrompt_mu_iso0p1,*den_medPrompt_mu);
  TEfficiency *eff_medPrompt_mu_iso0p2  = new TEfficiency(*num_SUSY_medPrompt_mu_iso0p2,*den_medPrompt_mu);
  TEfficiency *eff_medPrompt_mu_iso4  = new TEfficiency(*num_SUSY_medPrompt_mu_iso4,*den_medPrompt_mu);
  TEfficiency *eff_medPrompt_mu_iso6  = new TEfficiency(*num_SUSY_medPrompt_mu_iso6,*den_medPrompt_mu);

  //eff_medPrompt_mu_iso.push_back(eff_mp_mu); label_medPrompt_mu_iso.push_back("medium prompt (no iso)");
  eff_medPrompt_mu_iso.push_back(eff_medPrompt_mu_iso0p1); label_medPrompt_mu_iso.push_back("miniIso < 0.1");
  eff_medPrompt_mu_iso.push_back(eff_medPrompt_mu_iso0p2); label_medPrompt_mu_iso.push_back("miniIso < 0.2");
  eff_medPrompt_mu_iso.push_back(eff_medPrompt_mu_iso4); label_medPrompt_mu_iso.push_back("miniIso #upoint p_{T} < 4");
  eff_medPrompt_mu_iso.push_back(eff_medPrompt_mu_iso6); label_medPrompt_mu_iso.push_back("miniIso #upoint p_{T} < 6");

  PlotEff(eff_medPrompt_mu_iso,label_medPrompt_mu_iso,"eff_medPrompt_mu_iso", mu_label, ylabel_eff, samplename, true);

  //Muon medium prompt WP isolation fakerate
  vector<TEfficiency*> fr_medPrompt_mu_iso;

  TEfficiency *fr_medPrompt_mu_iso0p1  = new TEfficiency(*num_SUSY_medPrompt_mu_iso0p1_fr,*den_medPrompt_fake_mu);
  TEfficiency *fr_medPrompt_mu_iso0p2  = new TEfficiency(*num_SUSY_medPrompt_mu_iso0p2_fr,*den_medPrompt_fake_mu);
  TEfficiency *fr_medPrompt_mu_iso4  = new TEfficiency(*num_SUSY_medPrompt_mu_iso4_fr,*den_medPrompt_fake_mu);
  TEfficiency *fr_medPrompt_mu_iso6  = new TEfficiency(*num_SUSY_medPrompt_mu_iso6_fr,*den_medPrompt_fake_mu);

  //fr_medPrompt_mu_iso.push_back(fr_mp_mu); 
  fr_medPrompt_mu_iso.push_back(fr_medPrompt_mu_iso0p1); 
  fr_medPrompt_mu_iso.push_back(fr_medPrompt_mu_iso0p2); 
  fr_medPrompt_mu_iso.push_back(fr_medPrompt_mu_iso4); 
  fr_medPrompt_mu_iso.push_back(fr_medPrompt_mu_iso6); 

  PlotEff(fr_medPrompt_mu_iso,label_medPrompt_mu_iso,"fr_medPrompt_mu_iso", mu_label, ylabel_fr, samplename, false);

  //Other Plots
  /*
  Plot1D(dR_ele,"deltaR_ele",samplename);
  Plot1D(dR_mu, "deltaR_mu",samplename);

  Plot2D(miniIsoVsPt_sig,"miniIsoVsPt_sig_ele",ele_label,"MiniIso",samplename);
  Plot2D(miniIsoPtVsPt_sig,"miniIsoPtVsPt_sig_ele",ele_label,"MiniIso #upoint p_{T}",samplename);
  Plot2D(miniIsoVsPt_fake,"miniIsoVsPt_fake_ele",ele_label,"MiniIso",samplename);
  Plot2D(miniIsoPtVsPt_fake,"miniIsoPtVsPt_fake_ele",ele_label,"MiniIso #upoint p_{T}",samplename);

  Plot2D(miniIsoVsPt_mu_sig,"miniIsoVsPt_sig_mu",mu_label,"MiniIso",samplename);
  Plot2D(miniIsoPtVsPt_mu_sig,"miniIsoPtVsPt_sig_mu",mu_label,"MiniIso #upoint p_{T}",samplename);
  Plot2D(miniIsoVsPt_mu_fake,"miniIsoVsPt_fake_mu",mu_label,"MiniIso",samplename);
  Plot2D(miniIsoPtVsPt_mu_fake,"miniIsoPtVsPt_fake_mu",mu_label,"MiniIso #upoint p_{T}",samplename);
  */
  RocCurve(eff_sel, fake_sel, label_sel, nbins, "ROC_sel",samplename);
  RocCurve(eff_wp_ele, fr_wp_ele, label_wp_ele, nbins, "ROC_SUSY_wp_ele", samplename);
  RocCurve(eff_wp_v1, fr_wp_v1, label_wp_egamma, nbins, "ROC_Egamma_wp_v1", samplename);
  RocCurve(eff_wp_v2, fr_wp_v2, label_wp_egamma, nbins, "ROC_Egamma_wp_v2", samplename);
  RocCurve(eff_wp_mu, fr_wp_mu, label_wp_mu, nbins, "ROC_wp_mu", samplename);
  RocCurve(eff_sel_mu, fake_sel_mu, label_sel_mu, nbins, "ROC_sel_mu", samplename);

  return 0;
}

vector<int> clrs = {kBlue+2,kGreen+2,kRed+2,kOrange-3,kMagenta-1,kBlack,kAzure+2,kGray};
vector<int> markers = {20,21,22,23,29,33,34,43,49};

void PlotEff(vector<TEfficiency*> vec_eff, vector<TString> vec_names, TString name, TString xlabel, TString ylabel, TString sample, bool leg_pos, float ymax, float ymin){

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TMultiGraph* mg = new TMultiGraph();
  TCanvas *can = new TCanvas("can_"+name,"can_"+name,800,600);
  TLegend* leg;
  can->SetLeftMargin(0.15);
  can->SetGrid();

  for(int i = 0; i < vec_eff.size(); i++){
    vec_eff[i]->SetMarkerStyle(markers[i]);
    vec_eff[i]->SetMarkerColor(clrs[i]);
    vec_eff[i]->SetLineColor(clrs[i]);
    vec_eff[i]->SetStatisticOption(TEfficiency::kFNormal);
    mg->Add(vec_eff[i]->CreateGraph());
  }

  //vec_eff[0]->SetTitle(axis);
  if(leg_pos){
    leg = new TLegend(0.7,0.18,0.89,0.32);
  }
  if(!leg_pos) {
    leg = new TLegend(0.18,0.75,0.42,0.89);
  }

  mg->Draw("ap");
  mg->SetMinimum(ymin); 
  mg->SetMaximum(ymax); 
  mg->GetXaxis()->SetTitle(xlabel);
  mg->GetXaxis()->SetTitleOffset(1.25);
  mg->GetYaxis()->SetTitle(ylabel);
  leg->AddEntry(vec_eff[0],vec_names[0],"lep");

  for(int i = 1; i < vec_eff.size(); i++){
    vec_eff[i]->Draw("same");
    leg->AddEntry(vec_eff[i],vec_names[i],"lep");
    if(vec_eff.size() > 4)
      leg->SetNColumns(2);
    
  } 

  leg->Draw("same");

  CMSmark(RmSpace(sample));

  can->SaveAs("plots/"+name+"_"+sample+".pdf");
  can->Close();

  delete can;
  delete leg;
}

void Plot1D(TH1F* h1, TString name, TString sample){

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TCanvas *can = new TCanvas("can_"+name,"can_"+name,800,600);
  can->SetLogy();

  h1->GetXaxis()->SetTitle("#Delta R");
  h1->GetYaxis()->SetTitle("Events");

  h1->Draw();

  CMSmark(RmSpace(sample));

  can->SaveAs("plots/"+name+"_"+sample+".pdf");
}

void Plot2D(TH2F* h2, TString name, TString xlabel, TString ylabel, TString sample){

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TCanvas *can = new TCanvas("can_"+name,"can_"+name,800,600);
  can->SetLogz();

  h2->GetXaxis()->SetTitle(xlabel);
  h2->GetYaxis()->SetTitle(ylabel);

  h2->Draw("colz");

  CMSmark(RmSpace(sample));

  can->SaveAs("plots/"+name+"_"+sample+".pdf");
}

void RocCurve(vector<TEfficiency*> eff, vector<TEfficiency*> fakes, vector<TString> vec_names, int nbins, TString name, TString sample){

  if(eff.size() != fakes.size()){
    cout << "RocCurve: Efficiency vector must be same length as vector of fakes!" << endl;
    return;
  }

  int const size = eff.size();

  //vector<int> clrs = {kBlue+2,kGreen+2,kRed+2,kOrange-3,kMagenta-2,kBlack,kAzure+2,kGray};

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TCanvas *can = new TCanvas("can_"+name,"can_"+name,800,600);
  //TLegend* leg = new TLegend(0.18,0.75,0.32,0.89);
  //TLegend *leg = new TLegend(0.75,0.18,0.89,0.32);
  TLegend* leg = new TLegend(0.18,0.18,0.32,0.32); 

  can->SetLeftMargin(0.15);
  can->SetGrid();
  can->Draw();
  can->cd();

  vector<TGraphAsymmErrors*> tgae_vec;

  TMultiGraph* mg = new TMultiGraph();

  for (int i = 0; i < size; i++){
    tgae_vec.push_back(RocFromEff(eff[i], fakes[i], nbins));
    tgae_vec[i]->SetMarkerStyle(markers[i]);
    tgae_vec[i]->SetMarkerColor(clrs[i]);
    tgae_vec[i]->SetLineColor(clrs[i]);
    tgae_vec[i]->SetLineWidth(1);
    tgae_vec[i]->SetLineStyle(6);
    mg->Add(tgae_vec[i]);
  }

  mg->Draw("ap");
  mg->GetYaxis()->SetTitle("1-#varepsilon_{fake}");
  mg->GetXaxis()->SetTitle("#varepsilon_{true}");

  for (int i = 0; i < size; i++){
    leg->AddEntry(tgae_vec[i],vec_names[i],"lep");
  }

  leg->Draw("same");

  CMSmark(RmSpace(sample));
  can->SaveAs("plots/"+name+"_"+sample+".pdf");

}

TGraphAsymmErrors* RocFromEff(TEfficiency* eff, TEfficiency* fakes, int nbins){

  float *xAxis   = new float[nbins+1];
  float *yAxis   = new float[nbins+1]; 
  float *xErrLow = new float[nbins+1];
  float *xErrUp  = new float[nbins+1];
  float *yErrLow = new float[nbins+1];
  float *yErrUp  = new float[nbins+1];

  int N = 0;

  for(int i = 0; i < nbins; i++){
    
    if(fakes->GetEfficiency(i+1) <= 0. || eff->GetEfficiency(i+1) <= 0. ||
       fakes->GetEfficiency(i+1) > 1 || eff->GetEfficiency(i+1) > 1)
      continue;

    yAxis[N]   = 1-fakes->GetEfficiency(i+1);
    yErrLow[N] = fakes->GetEfficiencyErrorLow(i+1);
    yErrUp[N]  = fakes->GetEfficiencyErrorUp(i+1);

    xAxis[N]   = eff->GetEfficiency(i+1);
    xErrLow[N] = eff->GetEfficiencyErrorLow(i+1);
    xErrUp[N]  = eff->GetEfficiencyErrorUp(i+1);
  
    N++;

    // cout << i << " | " << nbins << " " << endl;
    // cout << xAxis[N-1] << " " << xErrLow[N-1] << " " << xErrUp[N-1] << endl;
    // cout << yAxis[N-1] << " " << yErrLow[N-1] << " " << yErrUp[N-1] << endl;
  }

  TGraphAsymmErrors *tgae = new TGraphAsymmErrors(N, xAxis, yAxis, xErrLow, xErrUp, yErrLow, yErrUp);

  delete [] xAxis;
  delete [] yAxis;
  delete [] xErrLow;
  delete [] xErrUp;
  delete [] yErrLow;
  delete [] yErrUp;

  return tgae;
}

void CMSmark(TString plotTitle){
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();
  l.SetTextSize(0.035);
  l.SetTextFont(42);
  l.DrawLatex(1 - (gPad->GetRightMargin() + plotTitle.Sizeof()*0.012), 1 - (gPad->GetTopMargin() - 0.017),plotTitle);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(0.15,0.91,"#bf{CMS} #it{Preliminary}");
}

TString RmSpace(TString str){

  TString temp_str = str;
  temp_str.ReplaceAll("_", " ");
  temp_str.ReplaceAll(".", "p");

  return temp_str;
}

