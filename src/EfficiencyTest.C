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

enum lepType {kSignal,kTau,kConversion,kLight,kHeavyB,kHeavyC,kUnmatched};

void PlotEff(vector<TEfficiency*> vec_eff, vector<TString> vec_names, TString name, TString xlabel, TString ylabel, TString sample, bool isEff, float ymax = 1.05, float ymin = 0.);
void Plot1D(TH1F* h1, TString name, TString sample);
void Plot2D(TH2F* h2, TString name, TString xlabel, TString ylabel, TString sample);
void RocCurve(vector<TEfficiency*> eff, vector<TEfficiency*> fakes, vector<TString> vec_names, int nbins, TString name, TString sample);
void PlotCategories(TH1F* denominator, vector<TH1F*> numerators, vector<TString> label_lep, TString name,
                    TString xlabel, TString ylabel, TString samplename, bool leg_pos, float xmax=1.05, float xmin=0.);
void CMSmark(TString plotTitle);

TString RmSpace(TString str);
TGraphAsymmErrors* RocFromEff(TEfficiency* eff, TEfficiency* fakes, int nbins);

void InitHist(vector<TH1F*> &histArr, int nArr, int nBins, float xmin, float xmax, TString name);
int GetMatchedMomID(TLorentzVector leptonCand, bool &isMatched, bool &isSignal, KUAnalysis *ku);
int AssignLeptonMomType(int motherID, bool isMatched, bool isSignal);
TH1F *AddHists(vector<TH1F*> histArr, TString name);

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

  int nbins = 40;
  float xmin = 0.;
  float xmax = 100.;

  TString ylabel_eff = "Efficiency";
  TString ylabel_fr = "Fakerate";
  TString ylabel2 = "Fractional Contribution";
  TString ele_label = "Electron p_{T} [GeV]";
  TString mu_label = "Muon p_{T} [GeV]";

  int nCat = 7;
  //electrons
  vector<TH1F*> den_ele_cat; InitHist(den_ele_cat, nCat, nbins, xmin, xmax, "den_ele_cat");
  vector<TH1F*> num_ele_bronze; InitHist(num_ele_bronze, nCat, nbins, xmin, xmax, "num_ele_bronze");
  vector<TH1F*> num_ele_silver; InitHist(num_ele_silver, nCat, nbins, xmin, xmax, "num_ele_silver");
  vector<TH1F*> num_ele_gold; InitHist(num_ele_gold, nCat, nbins, xmin, xmax, "num_ele_gold");
  vector<TH1F*> num_ele_sip3DnoIso; InitHist(num_ele_sip3DnoIso, nCat, nbins, xmin, xmax, "num_ele_sip3DnoIso");
  vector<TH1F*> num_ele_isoNoSip3D; InitHist(num_ele_isoNoSip3D, nCat, nbins, xmin, xmax, "num_ele_isoNoSip3D");
  vector<TH1F*> num_ele_noIsoNoSip; InitHist(num_ele_noIsoNoSip, nCat, nbins, xmin, xmax, "num_ele_noIsoNoSip");
  //muons
  vector<TH1F*> den_mu_cat; InitHist(den_mu_cat, nCat, nbins, xmin, xmax, "den_mu_cat");
  vector<TH1F*> num_mu_bronze; InitHist(num_mu_bronze, nCat, nbins, xmin, xmax, "num_mu_bronze");
  vector<TH1F*> num_mu_silver; InitHist(num_mu_silver, nCat, nbins, xmin, xmax, "num_mu_silver");
  vector<TH1F*> num_mu_gold; InitHist(num_mu_gold, nCat, nbins, xmin, xmax, "num_mu_gold");
  vector<TH1F*> num_mu_sip3DnoIso; InitHist(num_mu_sip3DnoIso, nCat, nbins, xmin, xmax, "num_mu_sip3DnoIso");
  vector<TH1F*> num_mu_isoNoSip3D; InitHist(num_mu_isoNoSip3D, nCat, nbins, xmin, xmax, "num_mu_isoNoSip3D");
  vector<TH1F*> num_mu_noIsoNoSip; InitHist(num_mu_noIsoNoSip, nCat, nbins, xmin, xmax, "num_mu_noIsoNoSip");

  int nEleCand = -999;
  int nMuCand = -999;
  int eleMomID, muMomID;
  int momType = -999;

  float weight = -999.;
  float pt = -999.;
  float miniIso = -999.;
  float sip3D = -999.;
  float deltaR = -999.;

  bool isMatched = false;
  bool isSignal = false;

  TLorentzVector gen_part, reco_ele, reco_mu;
  vector<TLorentzVector> gen_ele_vec, gen_mu_vec;

  string line;
  TChain *tch = new TChain("KUAnalysis");

  ifstream file (txtfile);
  if (file.is_open()){
    while (!file.eof()){
      getline(file,line);
      tch->Add(TString(line));
    }
  }
  else
    tch->Add(filename);

  KUAnalysis *ku = new KUAnalysis(tch);

  int nEntries = ku->fChain->GetEntries();
  if(test)
    nEntries = 100000;

  //Loop over entries in root file                                                                                                                  
  for(int e = 0; e < nEntries; e++){
    if (e % 10000 == 0) {
      fprintf(stdout, "\r  Processed events: %8d of %8d ", e, nEntries);
    }
    fflush(stdout);

    ku->fChain->GetEntry(e);

    weight = 1.0;//ku->weight;
    nEleCand = ku->Nele;

    //cout << "Event: " << e << endl;
    for(int r = 0; r < ku->Nele; r++){
      pt = ku->Electron_pt->at(r);
      miniIso = ku->MiniIso_ele->at(r);
      sip3D = ku->Electron_sip3d->at(r);
      deltaR = 999.;

      reco_ele.SetPtEtaPhiM(ku->PT_ele->at(r),ku->Eta_ele->at(r),ku->Phi_ele->at(r),ku->M_ele->at(r));
      //cout << " electron candidate: " << r << endl; 
      eleMomID = GetMatchedMomID(reco_ele, isMatched, isSignal, ku);
      momType = AssignLeptonMomType(eleMomID, isMatched, isSignal);
      if(ku->ID_ele->at(r) >= 1 && pt > 5){
	den_ele_cat[momType]->Fill(pt,weight);
	//bronze electrons                                                                                                               
	if(ku->ID_ele->at(r) >= 1 && ku->ID_ele->at(r) < 4 && sip3D >= 4 && miniIso*pt >= 5)
	  num_ele_bronze[momType]->Fill(pt,weight);
	//silver electrons                                                                                                               
	if(ku->ID_ele->at(r) >= 4 && sip3D >= 4 && miniIso*pt < 5)
	  num_ele_silver[momType]->Fill(pt,weight);
	//gold electrons                                                                                                                 
	if(ku->ID_ele->at(r) >= 4 && sip3D < 4 && miniIso*pt < 5)
	  num_ele_gold[momType]->Fill(pt,weight);
	//other electron categories                                                                                                               
	if(ku->ID_ele->at(r) >= 4 && sip3D < 4)
	  num_ele_sip3DnoIso[momType]->Fill(pt,weight);
	if(ku->ID_ele->at(r) >= 4 && miniIso*pt < 5)
	  num_ele_isoNoSip3D[momType]->Fill(pt,weight);
	if(ku->ID_ele->at(r) >= 4)
	  num_ele_noIsoNoSip[momType]->Fill(pt,weight);
      }
    }

    for(int r = 0; r < ku->Nmu; r++){
      pt = ku->Muon_pt->at(r);
      miniIso = ku->MiniIso_mu->at(r);
      sip3D = ku->Muon_sip3d->at(r);

      reco_mu.SetPtEtaPhiM(ku->PT_mu->at(r),ku->Eta_mu->at(r),ku->Phi_mu->at(r),ku->M_mu->at(r));
      //cout << " muon candidate: " << r << endl;
      muMomID = GetMatchedMomID(reco_mu, isMatched, isSignal, ku);
      momType = AssignLeptonMomType(muMomID, isMatched, isSignal);

      if(ku->ID_mu->at(r) >= 0 && pt > 5){
        den_mu_cat[momType]->Fill(pt,weight);
        //bronze muons
        if(ku->Muon_mediumId->at(r) && ku->ID_mu->at(r) < 4 && sip3D >= 4 && miniIso*pt >= 5)
          num_mu_bronze[momType]->Fill(pt,weight);
        //silver muons
        if(ku->Muon_mediumId->at(r) && sip3D >= 4 && miniIso*pt < 5)
          num_mu_silver[momType]->Fill(pt,weight);
        //gold muons
        if(ku->Muon_mediumId->at(r) && sip3D < 4 && miniIso*pt < 5)
          num_mu_gold[momType]->Fill(pt,weight);
        //other muons
        if(ku->Muon_mediumId->at(r) && sip3D < 4)
          num_mu_sip3DnoIso[momType]->Fill(pt,weight);
        if(ku->Muon_mediumId->at(r) && miniIso*pt < 5)
          num_mu_isoNoSip3D[momType]->Fill(pt,weight);
        if(ku->Muon_mediumId->at(r))
          num_mu_noIsoNoSip[momType]->Fill(pt,weight);
      }
    }

  }

  vector<TString> leg_label;
  leg_label.push_back("Signal");
  leg_label.push_back("#tau lepton");
  leg_label.push_back("Conversion");
  leg_label.push_back("Light");
  leg_label.push_back("Heavy b");
  leg_label.push_back("Heavy c");
  leg_label.push_back("Unmatched");

  //electron 
  TH1F *den_ele_bronze     = AddHists(num_ele_bronze,"den_ele_bronze");
  TH1F *den_ele_silver     = AddHists(num_ele_silver,"den_ele_silver");
  TH1F *den_ele_gold       = AddHists(num_ele_gold,"den_ele_gold");
  TH1F *den_ele_sip3DnoIso = AddHists(num_ele_sip3DnoIso,"den_ele_sip3DnoIso");
  TH1F *den_ele_isoNoSip3D = AddHists(num_ele_isoNoSip3D,"den_ele_isoNoSip3D");
  TH1F *den_ele_noIsoNoSip = AddHists(num_ele_noIsoNoSip,"den_ele_noIsoNoSip");

  PlotCategories(den_ele_bronze, num_ele_bronze, leg_label, "bronze_ele", ele_label, ylabel2, samplename, false);
  PlotCategories(den_ele_silver, num_ele_silver, leg_label, "silver_ele", ele_label, ylabel2, samplename, false);
  PlotCategories(den_ele_gold, num_ele_gold, leg_label, "gold_ele", ele_label, ylabel2, samplename, true);
  PlotCategories(den_ele_sip3DnoIso, num_ele_sip3DnoIso, leg_label, "sip3DnoIso_ele", ele_label, ylabel2, samplename, true);
  PlotCategories(den_ele_isoNoSip3D, num_ele_isoNoSip3D, leg_label, "isoNoSip3D_ele", ele_label, ylabel2, samplename, true);
  PlotCategories(den_ele_noIsoNoSip, num_ele_noIsoNoSip, leg_label, "noIsoNoSip_ele", ele_label, ylabel2, samplename, false);

  //muon
  TH1F *den_mu_bronze      = AddHists(num_mu_bronze,"den_mu_bronze");
  TH1F *den_mu_silver      = AddHists(num_mu_silver,"den_mu_silver");
  TH1F *den_mu_gold        = AddHists(num_mu_gold,"den_mu_gold");
  TH1F *den_mu_sip3DnoIso  = AddHists(num_mu_sip3DnoIso,"den_mu_sip3DnoIso");
  TH1F *den_mu_isoNoSip3D  = AddHists(num_mu_isoNoSip3D,"den_mu_isoNoSip3D");
  TH1F *den_mu_noIsoNoSip  = AddHists(num_mu_noIsoNoSip,"den_mu_noIsoNoSip");

  PlotCategories(den_mu_bronze, num_mu_bronze, leg_label, "bronze_mu", mu_label, ylabel2, samplename, false);
  PlotCategories(den_mu_silver, num_mu_silver, leg_label, "silver_mu", mu_label, ylabel2, samplename, false);
  PlotCategories(den_mu_gold, num_mu_gold, leg_label, "gold_mu", mu_label, ylabel2, samplename, true);
  PlotCategories(den_mu_sip3DnoIso, num_mu_sip3DnoIso, leg_label, "sip3DnoIso_mu", mu_label, ylabel2, samplename, true);
  PlotCategories(den_mu_isoNoSip3D, num_mu_isoNoSip3D, leg_label, "isoNoSip3D_mu", mu_label, ylabel2, samplename, true);
  PlotCategories(den_mu_noIsoNoSip, num_mu_noIsoNoSip, leg_label, "noIsoNoSip_mu", mu_label, ylabel2, samplename, false);

  return 0;
}

vector<int> clrs = {kBlue+2,kGreen+2,kRed+2,kOrange-3,kMagenta-1,kBlack,kAzure+2,kGray};
vector<int> markers = {20,21,22,23,29,33,34,43,49};

void PlotEff(vector<TEfficiency*> vec_eff, vector<TString> vec_names, TString name, 
	     TString xlabel, TString ylabel, TString sample, bool leg_pos, float ymax, float ymin){

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

  //can->SetLogx();
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

void PlotCategories(TH1F* denominator, vector<TH1F*> numerators, vector<TString> label_lep, TString name, 
		    TString xlabel, TString ylabel, TString samplename, bool leg_pos, float xmax, float xmin){

  int nNum = numerators.size();

  vector<TEfficiency*> eff_lep;

  for(int i = 0; i < nNum; i++){
    eff_lep.push_back(new TEfficiency(*numerators[i],*denominator));
  }

  PlotEff(eff_lep, label_lep, name, xlabel, ylabel, samplename, leg_pos);

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
  l.DrawLatex(0.15,0.91,"#bf{CMS} #it{Simulation}");
}

TString RmSpace(TString str){

  TString temp_str = str;
  temp_str.ReplaceAll("_", " ");
  temp_str.ReplaceAll(".", "p");

  return temp_str;
}

void InitHist(vector<TH1F*> &histArr, int nArr, int nBins, float xmin, float xmax, TString name){
  TString temp_name;

  for(int i = 0; i < nArr; i++){
    temp_name = name+to_string(i);
    histArr.push_back(new TH1F(temp_name, temp_name, nBins, xmin, xmax));
  }
}

int GetMatchedMomID(TLorentzVector leptonCand, bool &isMatched, bool &isSignal, KUAnalysis *ku){

  int genPartID = -999;
  int genPartMotherID = -999;
  int genPartMotherIdx = -999;
  int leptonMomID = -999;
  int leptonStatus = -999;
  int momStatus = -999;

  float deltaR = 999.;

  TLorentzVector gen_part;

  isMatched = false;
  isSignal = false;

  for(int g = 0; g < ku->nGenPart; g++){
    gen_part.SetPtEtaPhiM(ku->GenPart_pt->at(g),ku->GenPart_eta->at(g),ku->GenPart_phi->at(g),ku->GenPart_mass->at(g));

    //if(ku->GenPart_status->at(g) == 23 && ku->GenPart_pdgId->at(ku->GenPart_genPartIdxMother->at(g)))
    //cout << ku->GenPart_pdgId->at(ku->GenPart_genPartIdxMother->at(g)) << endl;
    leptonStatus = ku->GenPart_status->at(g);
    genPartID = ku->GenPart_pdgId->at(g);
    if(leptonCand.DeltaR(gen_part) >= deltaR || (abs(genPartID) != 11 && abs(genPartID) != 13) || leptonStatus != 1)
      continue;

    deltaR = leptonCand.DeltaR(gen_part);

    genPartMotherIdx = ku->GenPart_genPartIdxMother->at(g);
    if(genPartMotherIdx >= 0 && deltaR < 0.05){
      //cout << "   DeltaR: " << deltaR << endl;                                                                             
      //cout << "   gen ID: " << genPartID << endl;
      genPartMotherID = ku->GenPart_pdgId->at(genPartMotherIdx);
      while(abs(genPartMotherID) == 11 || abs(genPartMotherID) == 13 && genPartMotherIdx >= 0){
	momStatus = ku->GenPart_status->at(genPartMotherIdx);
	if(momStatus == 23 && (abs(genPartMotherID) == 11 || abs(genPartMotherID) == 13)){
	   isSignal = true;
	   isMatched = true;
	   return genPartMotherID;
	}
	//cout << genPartMotherID << endl;
	genPartMotherIdx = ku->GenPart_genPartIdxMother->at(genPartMotherIdx);
	if(genPartMotherIdx >= 0)
	  genPartMotherID = ku->GenPart_pdgId->at(genPartMotherIdx);
	else continue;
      }
      cout << "   gen Mom ID: " << genPartMotherID << endl;
      if(abs(genPartMotherID) != 11 && abs(genPartMotherID) != 13 && !isSignal){
	isMatched = true;
	leptonMomID = genPartMotherID;
	//cout << "   Matched to electron " << r << endl;
      }
    }
  }
  return leptonMomID;
}

  int AssignLeptonMomType(int motherID, bool isMatched, bool isSignal){

  lepType type;

  //cout << motherID << endl;
  if(isMatched){
    //if(abs(motherID) == 24 || motherID == 23)
    if(isSignal || abs(motherID) == 24 || abs(motherID) == 23)
      type = kSignal;
    else if(abs(motherID) == 15)
      type = kTau;
    else if((abs(motherID%1000) > 100 && abs(motherID%1000) < 400)
       || (abs(motherID%1000) > 1000 && abs(motherID%1000) < 4000)
       || (abs(motherID) > 0 && abs(motherID) < 4)
       || motherID == 21)
      type = kLight;
    else if((abs(motherID%1000) > 400 && abs(motherID%1000) < 500)
       || (abs(motherID%1000) > 4000 && abs(motherID%1000) < 5000)
       || (abs(motherID) > 3 && abs(motherID) < 5))
      type = kHeavyC;
    else if((abs(motherID%1000) > 500 && abs(motherID%1000) < 600)
       || (abs(motherID%1000) > 5000 && abs(motherID%1000) < 6000)
       || (abs(motherID) > 4 && abs(motherID) < 7))
      type = kHeavyB;
    else //(motherID == 22)
      type = kConversion;

  }
  else
    type = kUnmatched;
  return type;
}

TH1F *AddHists(vector<TH1F*> histArr, TString name){

  int nHist = histArr.size();
  TH1F *combinedHist = (TH1F*)histArr[0]->Clone(name);

  for(int i = 1; i < nHist; i++)
    combinedHist->Add(histArr[i]);
  
  return combinedHist;
}
