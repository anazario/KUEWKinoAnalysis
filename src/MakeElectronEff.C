#include <iostream>
#include <vector>

#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLorentzVector.h>

#include "KUAnalysis.h"

using namespace std;

void PlotEff(vector<TEfficiency*> vec_eff, vector<TString> vec_names, TString name, bool isEff);
void Plot1D(TH1F* h1);

int main(int argc, char *argv[]){

  TLorentzVector tmp_gen_ele, ele_cand;
  vector<TLorentzVector> gen_ele;
  vector<int> passIdx;
  int idx = 0;

  int nbins = 40;
  float xmin = 0.;
  float xmax = 100.;
  float dPhi_cut = 0.04;
  float ip3D = 0.0175;
  float sip3D = 2.5;
  float dxy = 0.05;
  float dz = 0.1;

  bool isGenZ;

  TH1F *dPhi = new TH1F("dPhi","dPhi", 100, 0., 1.);

  //Working points (SUSY)
  TH1F *num_vl = new TH1F("num_vl","num_vl",nbins,xmin,xmax);
  TH1F *num_l = new TH1F("num_l","num_l",nbins,xmin,xmax);
  TH1F *num_m = new TH1F("num_m","num_m",nbins,xmin,xmax);
  TH1F *num_t = new TH1F("num_t","num_t",nbins,xmin,xmax);
  
  //Working points V1 (EGamma)
  TH1F *num_WP80V1 = new TH1F("num_WP80V1","num_WP80V1",nbins,xmin,xmax);
  TH1F *num_WP90V1 = new TH1F("num_WP90V1","num_WP90V1",nbins,xmin,xmax);
  TH1F *num_WPLV1 = new TH1F("num_WPLV1","num_WPLV1",nbins,xmin,xmax);

  //Working points V2 (EGamma)
  TH1F *num_WP80V2 = new TH1F("num_WP80V2","num_WP80V2",nbins,xmin,xmax);
  TH1F *num_WP90V2 = new TH1F("num_WP90V2","num_WP90V2",nbins,xmin,xmax);
  TH1F *num_WPLV2 = new TH1F("num_WPLV2","num_WPLV2",nbins,xmin,xmax);

  //Working points Spring16 (EGamma)
  TH1F *num_WP80S16 = new TH1F("num_WP80S16","num_WP80S16",nbins,xmin,xmax);
  TH1F *num_WP90S16 = new TH1F("num_WP90S16","num_WP90S16",nbins,xmin,xmax);
  TH1F *num_WPLS16 = new TH1F("num_WPLS16","num_WPLS16",nbins,xmin,xmax);

  //Electron candidate selection parameter histograms
  TH1F *num_lostHits = new TH1F("num_lostHits","num_lostHits",nbins,xmin,xmax);
  TH1F *num_ip3D = new TH1F("num_ip3D","num_ip3D",nbins,xmin,xmax);
  TH1F *num_sip3D = new TH1F("num_sip3D","num_sip3D",nbins,xmin,xmax);
  TH1F *num_convVeto = new TH1F("num_convVeto","num_convVeto",nbins,xmin,xmax);
  TH1F *num_dxy = new TH1F("num_dxy","num_dxy",nbins,xmin,xmax);
  TH1F *num_dz = new TH1F("num_dz","num_dz",nbins,xmin,xmax);

  //Fakerate
  TH1F *num_lostHits_f = new TH1F("num_lostHits_f","num_lostHits_f",nbins,xmin,xmax);
  TH1F *num_ip3D_f = new TH1F("num_ip3D_f","num_ip3D_f",nbins,xmin,xmax);
  TH1F *num_sip3D_f = new TH1F("num_sip3D_f","num_sip3D_f",nbins,xmin,xmax);
  TH1F *num_convVeto_f = new TH1F("num_convVeto_f","num_convVeto_f",nbins,xmin,xmax);
  TH1F *num_dxy_f = new TH1F("num_dxy_f","num_dxy_f",nbins,xmin,xmax);
  TH1F *num_dz_f = new TH1F("num_dz_f","num_dz_f",nbins,xmin,xmax);

  //Denominator (minimum electron candidate requirements)
  TH1F *den_all = new TH1F("den_all","den_all",nbins,xmin,xmax);
  TH1F *den_fake = new TH1F("den_fake","den_fake",nbins,xmin,xmax);

  TChain *tch = new TChain("KUAnalysis");
  tch->Add("root/test.root");

  KUAnalysis *ku = new KUAnalysis(tch);

  int nEntries = ku->fChain->GetEntries();

  for(int e = 0; e < nEntries; e++){
    if (e % 1000 == 0) {
      fprintf(stdout, "\r  Processed events: %8d of %8d ", e, nEntries);
    }
    fflush(stdout);

    ku->fChain->GetEntry(e);

    isGenZ = false;

    for(int g = 0; g < ku->genNele; g++){
      if(ku->genMomPDGID_ele->at(g) == 23){
	isGenZ = true;
	tmp_gen_ele.SetPtEtaPhiM(ku->genPT_ele->at(g),ku->genEta_ele->at(g),ku->genPhi_ele->at(g),ku->genM_ele->at(g));
	gen_ele.push_back(tmp_gen_ele);
      }
    }

    if(isGenZ){

      for(int i = 0; i < ku->Nele; i++){

	ele_cand.SetPtEtaPhiM(ku->PT_ele->at(i),ku->Eta_ele->at(i),ku->Phi_ele->at(i),ku->M_ele->at(i));
	for(int g = 0; g < gen_ele.size(); g++){
	  dPhi->Fill(ele_cand.DeltaR(gen_ele[g]));
	  if(ele_cand.DeltaR(gen_ele[g]) < dPhi_cut){
	    passIdx.push_back(i);
	    break;
	  }
	}
      }
      gen_ele.clear();

      for(int i = 0; i < passIdx.size(); i++){ 

	idx = passIdx[i];

	if(ku->ID_ele->at(idx) >= 0){
	  den_all->Fill(ku->Electron_pt->at(idx));
	  //SUSY wp
	  if(ku->ID_ele->at(idx) >= 1)
	    num_vl->Fill(ku->Electron_pt->at(idx));
	  if(ku->ID_ele->at(idx) >= 2)
	    num_l->Fill(ku->Electron_pt->at(idx));
	  if(ku->ID_ele->at(idx) >= 3)
	    num_m->Fill(ku->Electron_pt->at(idx));
	  if(ku->ID_ele->at(idx) >= 4)
	    num_t->Fill(ku->Electron_pt->at(idx));

	  //EGamma wp Fall17V1
	  if(ku->Electron_mvaFall17V1noIso_WP80->at(idx))
            num_WP80V1->Fill(ku->Electron_pt->at(idx));
          if(ku->Electron_mvaFall17V1noIso_WP90->at(idx))
            num_WP90V1->Fill(ku->Electron_pt->at(idx));
          if(ku->Electron_mvaFall17V1noIso_WPL->at(idx))
            num_WPLV1->Fill(ku->Electron_pt->at(idx));
	  //EGamma wp Fall17V2
          if(ku->Electron_mvaFall17V2noIso_WP80->at(idx))
            num_WP80V2->Fill(ku->Electron_pt->at(idx));
          if(ku->Electron_mvaFall17V2noIso_WP90->at(idx))
            num_WP90V2->Fill(ku->Electron_pt->at(idx));
          if(ku->Electron_mvaFall17V2noIso_WPL->at(idx))
            num_WPLV2->Fill(ku->Electron_pt->at(idx));
	  //EGamma wp Spring16
          if(ku->Electron_mvaSpring16GP_WP80->at(idx))
            num_WP80S16->Fill(ku->Electron_pt->at(idx));
          if(ku->Electron_mvaSpring16GP_WP90->at(idx))
            num_WP90S16->Fill(ku->Electron_pt->at(idx));
          if(ku->Electron_mvaSpring16HZZ_WPL->at(idx))
            num_WPLS16->Fill(ku->Electron_pt->at(idx));

	  //electron candidate requirements
	  if(ku->Electron_lostHits->at(idx) == 0)
            num_lostHits->Fill(ku->Electron_pt->at(idx));
	  if(ku->Electron_ip3d->at(idx) < ip3D)
	    num_ip3D->Fill(ku->Electron_pt->at(idx));
	  if(ku->Electron_sip3d->at(idx) < sip3D)
	    num_sip3D->Fill(ku->Electron_pt->at(idx));
	  if(ku->Electron_convVeto->at(idx))
	    num_convVeto->Fill(ku->Electron_pt->at(idx));
	  if(fabs(ku->Electron_dxy->at(idx)) < dxy)
	    num_dxy->Fill(ku->Electron_pt->at(idx));
	  if(fabs(ku->Electron_dz->at(idx)) < dz)
	    num_dz->Fill(ku->Electron_pt->at(idx));
	}
      }
      passIdx.clear();      
    }
    
    if(!isGenZ){
      for(int i = 0; i < ku->Nele; i++){
	//Selection criteria fakerates
	if(ku->ID_ele->at(i) >= 0){
	  den_fake->Fill(ku->Electron_pt->at(i));
	  if(ku->Electron_lostHits->at(i) == 0)
	    num_lostHits_f->Fill(ku->Electron_pt->at(i));
	  if(ku->Electron_ip3d->at(i) < ip3D)
	    num_ip3D_f->Fill(ku->Electron_pt->at(i));
	  if(ku->Electron_sip3d->at(i) < sip3D)
	    num_sip3D_f->Fill(ku->Electron_pt->at(i));
	  if(ku->Electron_convVeto->at(i))
	    num_convVeto_f->Fill(ku->Electron_pt->at(i));
	  if(fabs(ku->Electron_dxy->at(i)) < dxy)
	    num_dxy_f->Fill(ku->Electron_pt->at(i));
	  if(fabs(ku->Electron_dz->at(i)) < dz)
	    num_dz_f->Fill(ku->Electron_pt->at(i));
	}
      }
    }
  }
  cout << endl;

  //Working point efficiencies (SUSY group)
  vector<TEfficiency*> eff_wp;
  vector<TString> label_wp;

  TEfficiency *eff_vl = new TEfficiency(*num_vl,*den_all);
  TEfficiency *eff_l = new TEfficiency(*num_l,*den_all);
  TEfficiency *eff_m = new TEfficiency(*num_m,*den_all);
  TEfficiency *eff_t = new TEfficiency(*num_t,*den_all);

  eff_wp.push_back(eff_vl); label_wp.push_back("baseline");
  eff_wp.push_back(eff_l); label_wp.push_back("very loose");
  eff_wp.push_back(eff_m); label_wp.push_back("loose");
  eff_wp.push_back(eff_t); label_wp.push_back("tight");

  PlotEff(eff_wp,label_wp,"eff_wp",true);

  //Working point efficiencies (Egamma Fall17V1)
  vector<TEfficiency*> eff_wp_v1;
  vector<TString> label_wp_egamma;

  TEfficiency *eff_wp80_v1 = new TEfficiency(*num_WP80V1,*den_all);
  TEfficiency *eff_wp90_v1 = new TEfficiency(*num_WP90V1,*den_all);
  TEfficiency *eff_wpl_v1 = new TEfficiency(*num_WPLV1,*den_all);

  eff_wp_v1.push_back(eff_wpl_v1); label_wp_egamma.push_back("WPL");  
  eff_wp_v1.push_back(eff_wp90_v1); label_wp_egamma.push_back("WP90");
  eff_wp_v1.push_back(eff_wp80_v1); label_wp_egamma.push_back("WP80");

  PlotEff(eff_wp_v1,label_wp_egamma,"eff_wp_egV1",true);

  //Working point efficiencies (Egamma Fall17V2)
  vector<TEfficiency*> eff_wp_v2;

  TEfficiency *eff_wp80_v2 = new TEfficiency(*num_WP80V2,*den_all);
  TEfficiency *eff_wp90_v2 = new TEfficiency(*num_WP90V2,*den_all);
  TEfficiency *eff_wpl_v2 = new TEfficiency(*num_WPLV2,*den_all);

  eff_wp_v2.push_back(eff_wpl_v2);
  eff_wp_v2.push_back(eff_wp90_v2); 
  eff_wp_v2.push_back(eff_wp80_v2);

  PlotEff(eff_wp_v2,label_wp_egamma,"eff_wp_egV2",true);

  //Working point efficiencies (Egamma Summer16)
  vector<TEfficiency*> eff_wp_s16;

  TEfficiency *eff_wp80_s16 = new TEfficiency(*num_WP80S16,*den_all);
  TEfficiency *eff_wp90_s16 = new TEfficiency(*num_WP90S16,*den_all);
  TEfficiency *eff_wpl_s16 = new TEfficiency(*num_WPLS16,*den_all);

  eff_wp_s16.push_back(eff_wpl_s16);
  eff_wp_s16.push_back(eff_wp90_s16);
  eff_wp_s16.push_back(eff_wp80_s16);

  PlotEff(eff_wp_s16,label_wp_egamma,"eff_wp_egS16",true);

  //electron candidate selection efficiencies and fakerates
  vector<TEfficiency*> eff_sel;
  vector<TEfficiency*> fake_sel;
  vector<TString> label_sel;

  TEfficiency *eff_lostHits = new TEfficiency(*num_lostHits,*den_all);
  TEfficiency *eff_ip3D = new TEfficiency(*num_ip3D,*den_all);
  TEfficiency *eff_sip3D = new TEfficiency(*num_sip3D,*den_all);
  TEfficiency *eff_convVeto = new TEfficiency(*num_convVeto,*den_all);
  TEfficiency *eff_dxy = new TEfficiency(*num_dxy,*den_all);
  TEfficiency *eff_dz = new TEfficiency(*num_dz,*den_all);

  eff_sel.push_back(eff_lostHits); label_sel.push_back(Form("pass lostHits"));
  eff_sel.push_back(eff_ip3D); label_sel.push_back(Form("ip3D < %.4f",ip3D));
  eff_sel.push_back(eff_sip3D); label_sel.push_back(Form("sip3D < %.1f",sip3D));
  eff_sel.push_back(eff_convVeto); label_sel.push_back("pass convVeto");
  eff_sel.push_back(eff_dxy); label_sel.push_back(Form("dxy < %.2f",dxy));
  eff_sel.push_back(eff_dz); label_sel.push_back(Form("dz < %.2f",dz));

  TEfficiency *fake_lostHits = new TEfficiency(*num_lostHits_f,*den_fake);
  TEfficiency *fake_ip3D = new TEfficiency(*num_ip3D_f,*den_fake);
  TEfficiency *fake_sip3D = new TEfficiency(*num_sip3D_f,*den_fake);
  TEfficiency *fake_convVeto = new TEfficiency(*num_convVeto_f,*den_fake);
  TEfficiency *fake_dxy = new TEfficiency(*num_dxy_f,*den_fake);
  TEfficiency *fake_dz = new TEfficiency(*num_dz_f,*den_fake);

  fake_sel.push_back(fake_lostHits);
  fake_sel.push_back(fake_ip3D);
  fake_sel.push_back(fake_sip3D);
  fake_sel.push_back(fake_convVeto);
  fake_sel.push_back(fake_dxy);
  fake_sel.push_back(fake_dz);

  PlotEff(eff_sel,label_sel,"eff_sel",true);
  PlotEff(fake_sel,label_sel,"fake_sel",false);

  Plot1D(dPhi);

  return 0;
}

void PlotEff(vector<TEfficiency*> vec_eff, vector<TString> vec_names, TString name, bool isEff){

  vector<int> clrs = {kBlue+2,kGreen+2,kRed+2,kOrange-3,kMagenta-2,kBlack,kAzure+2,kGray};

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TCanvas *can = new TCanvas("can_"+name,"can_"+name,800,600);
  TLegend* leg = new TLegend(0.75,0.18,0.89,0.32);
  can->SetLeftMargin(0.15);
  can->SetGrid();

  for(int i = 0; i < vec_eff.size(); i++){
    vec_eff[i]->SetMarkerStyle(10);
    vec_eff[i]->SetMarkerColor(clrs[i]);
    vec_eff[i]->SetLineColor(clrs[i]);
  }

  if(isEff)
    vec_eff[0]->SetTitle("; Electron p_{T} [GeV]; Efficiency");
  if(!isEff) 
    vec_eff[0]->SetTitle("; Electron p_{T} [GeV]; Fakerate");

  vec_eff[0]->Draw();
  leg->AddEntry(vec_eff[0],vec_names[0],"lep");

  for(int i = 1; i < vec_eff.size(); i++){
    vec_eff[i]->Draw("same");
    leg->AddEntry(vec_eff[i],vec_names[i],"lep");
  } 

  leg->Draw("same");

  gPad->Update(); 
  vec_eff[0]->GetPaintedGraph()->SetMinimum(0);
  gPad->Update(); 

  can->SaveAs("plots/"+name+".pdf");

}

void Plot1D(TH1F* h1){

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TCanvas *can = new TCanvas("can_","can_",800,600);
  can->SetLogy();

  h1->GetXaxis()->SetTitle("#Delta#phi");
  h1->GetYaxis()->SetTitle("Events");

  h1->Draw();

  can->SaveAs("plots/Dphi.pdf");
}
