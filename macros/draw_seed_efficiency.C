#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"

void draw_seed_efficiency(TString pid = "pi", TString dir1 = "../acts", TString dir2 = "../acts", TString dir3 = "../acts"){
  dir1.Append("/");
  dir2.Append("/");  
  dir3.Append("/");
  bool isPi = pid.Contains("pi");

  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadBottomMargin(0.11);

  TFile* f1 = new TFile(dir1 + "seed_efficiency.root");
  TFile* f2 = new TFile(dir2 + "seed_efficiency.root");
  TFile* f3 = new TFile(dir3 + "seed_efficiency.root");

  TH1D* hMcPt1 = (TH1D*) f1->Get(isPi ? "hMcPtPi16" : "hMcPtPr16");
  TH1D* hMcPt2 = (TH1D*) f2->Get(isPi ? "hMcPtPi19" : "hMcPtPr19");
  TH1D* hMcPt3 = (TH1D*) f3->Get(isPi ? "hMcPtPi16" : "hMcPtPr16");

  TH1D* hRcPt1 = (TH1D*) f1->Get(isPi ? "hRcPtPi16" : "hRcPtPr16");
  TH1D* hRcPt2 = (TH1D*) f2->Get(isPi ? "hRcPtPi19" : "hRcPtPr19");
  TH1D* hRcPt3 = (TH1D*) f3->Get(isPi ? "hRcPtPi16" : "hRcPtPr16");

  auto hEffPt1 = (TH1D*) hRcPt1->Clone("hEffPt1");
  auto hEffPt2 = (TH1D*) hRcPt2->Clone("hEffPt2");
  auto hEffPt3 = (TH1D*) hRcPt3->Clone("hEffPt3");

  hEffPt1->Divide(hRcPt1, hMcPt1, 1, 1, "B");
  hEffPt2->Divide(hRcPt2, hMcPt2, 1, 1, "B");
  hEffPt3->Divide(hRcPt3, hMcPt3, 1, 1, "B");  

  hEffPt1->SetTitle(";p_{T} (GeV);Efficiency");
  hEffPt1->SetLineWidth(2);
  hEffPt2->SetLineWidth(2);
  hEffPt3->SetLineWidth(2);
  hEffPt1->SetLineColor(kBlue);
  hEffPt2->SetLineColor(kMagenta);
  hEffPt3->SetLineColor(kRed);
  
  new TCanvas;
  hEffPt1->GetYaxis()->SetRangeUser(0,1.1);
  hEffPt1->SetLabelSize(0.045,"XY");
  hEffPt1->SetTitleSize(0.045,"XY");
  hEffPt1->SetTitleOffset(1.1,"X");
  hEffPt1->Draw();
  hEffPt2->Draw("same");
//  hEffPt3->Draw("same");


  TLine* l = new TLine();
  l->SetLineColor(kGray);
  l->DrawLine(0,1,1,1);

  TLegend* leg = new TLegend(0.7,0.40,0.9,0.6);
  leg->SetBorderSize(0);
  leg->AddEntry(hEffPt1,"#eta = 1.6");
  leg->AddEntry(hEffPt2,"#eta = 1.9");
//  leg->AddEntry(hEffPt3,"#eta = 2.2");
  leg->Draw();

  gPad->Print(dir1 + Form("seed_efficiency_%s.png",pid.Data()));
}
