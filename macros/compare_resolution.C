#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"

void compare_resolution(TString pid = "Pi"){
  bool isPi = pid.Contains("Pi");

  TFile* f1 = new TFile("single/resolution.root");
  TFile* f2 = new TFile("urqmd/resolution.root");

  TGraph* gRes1 = (TGraph*) f1->Get(Form("gRes%s1.60",pid.Data()));
  TGraph* gRes2 = (TGraph*) f2->Get(Form("gRes%s1.60",pid.Data()));

  TCanvas* c = new TCanvas("c", "c", 800, 600);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.07);
  TH1F* hFrame = gPad->DrawFrame(0.,0.,1.1,0.1);

  hFrame->SetTitle(";p_{T}^{MC} (GeV/c); p_{T} resolution");
  hFrame->SetTitle("#eta = 1.6");
  hFrame->GetXaxis()->SetTitleOffset(1.2);

  gRes1->SetLineStyle(9);
  gRes2->SetLineStyle(9);

  gRes1->SetLineWidth(3); gRes1->SetLineColor(kBlue);
  gRes2->SetLineWidth(3); gRes2->SetLineColor(kMagenta);

  gRes1->Draw("same");
  gRes2->Draw("same");

  gRes1->RemovePoint(0);
  gRes2->RemovePoint(0);

  TLegend* legend = new TLegend(0.35,0.68,0.70,0.85);
  legend->SetBorderSize(0);
  legend->AddEntry(gRes1,"Single tracks","l");
  legend->AddEntry(gRes2,"Central UrQMD","l");
  legend->Draw();
  
  gPad->Print(Form("resolution%s160.png",pid.Data()));
}

