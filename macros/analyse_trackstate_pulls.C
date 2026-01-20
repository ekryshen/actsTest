#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "vector"
#include "style.h"
#include "tree_trackstates.C"

void analyse_trackstate_pulls(){
  // TFile* f = new TFile("../build/trackstates_ckf.root");
  TFile* f = new TFile("../build/test/trackstates.root");
  TTree* t = (TTree*) f->Get("trackstates");
  SetBranchAddresses(t);

  TH1D* h_eQOP_fit = new TH1D("h_eQOP_fit","",200,0,2);
  TH1D* h_eLOC0_fit = new TH1D("h_eLOC0_fit","",200,-90,90);
  TH1D* h_eLOC1_fit = new TH1D("h_eLOC1_fit","",200,-90,90);
  TH1D* h_pull_eQOP_fit = new TH1D("h_pull_eQOP_fit","",200,-10,10);
  TH1D* h_pull_eLOC0_fit = new TH1D("h_pull_eLOC0_fit","",200,-10,10);
  TH1D* h_pull_eLOC1_fit = new TH1D("h_pull_eLOC1_fit","",200,-10,10);

  TH1D* h_pT_ubs = new TH1D("h_pT_ubs","",200,0,0.4);

  int nEvents = t->GetEntries();
  for (int ev=0; ev<nEvents; ev++){
    t->GetEntry(ev);
    // fill unbiased parameters and pulls at the first measurement
    int im = 0;
    if (m_hasParams[eUnbiased]->at(im)) {
      h_eQOP_fit->Fill(m_eQOP[eUnbiased]->at(im));
      h_eLOC0_fit->Fill(m_res_eLOC0[eUnbiased]->at(im));
      h_eLOC1_fit->Fill(m_res_eLOC1[eUnbiased]->at(im));
      h_pull_eQOP_fit->Fill(m_pull_eQOP[eUnbiased]->at(im));
      h_pull_eLOC0_fit->Fill(m_pull_eLOC0[eUnbiased]->at(im));
      h_pull_eLOC1_fit->Fill(m_pull_eLOC1[eUnbiased]->at(im));
      h_pT_ubs->Fill(m_pT[eUnbiased]->at(im));
    }
  }

  TCanvas* cPulls = new TCanvas("cPulls","",1900,800);
  cPulls->Divide(3,2,0.001,0.001);
  cPulls->cd(1);
  SetPad(gPad);
  SetHisto(h_eLOC0_fit,";d (mm)");
  h_eLOC0_fit->Draw();
  
  cPulls->cd(2);
  SetPad(gPad);
  SetHisto(h_eLOC1_fit,";z (mm)");
  h_eLOC1_fit->Draw();

  cPulls->cd(3);
  SetPad(gPad);
  SetHisto(h_eQOP_fit,";q/p (1/GeV)");
  h_eQOP_fit->Draw();
  
  cPulls->cd(4);
  SetPad(gPad);
  SetHisto(h_pull_eLOC0_fit,";pull d");
  h_pull_eLOC0_fit->Draw();

  cPulls->cd(5);
  SetPad(gPad);
  SetHisto(h_pull_eLOC1_fit,";pull z");
  h_pull_eLOC1_fit->Draw();

  cPulls->cd(6);
  SetPad(gPad);
  SetHisto(h_pull_eQOP_fit,";pull q/p");
  h_pull_eQOP_fit->Draw();
  h_pull_eQOP_fit->Fit("gaus","","",-2,2);

  cPulls->Print("pulls_mes0_unbiased.png");

  new TCanvas;
  h_pT_ubs->Draw();
}