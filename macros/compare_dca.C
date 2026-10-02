#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"

void compare_dca(TString pid = "Pi", double eta = 1.9){
  bool isPi = pid.Contains("Pi");

  TFile* f1 = new TFile("nazar2/merged/single/dca_resolution.root");
  TFile* f2 = new TFile("nazar2/merged/urqmd/dca_resolution.root");

  f1->ls();
  f2->ls();

  TGraph* gRes1 = (TGraph*) f1->Get(Form("gDcaRes%s%.2f",pid.Data(),eta));
  TGraph* gRes2 = (TGraph*) f2->Get(Form("gDcaRes%s%.2f",pid.Data(),eta));

  TCanvas* c = new TCanvas("c", "c", 800, 600);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.07);
  TH1F* hFrame = gPad->DrawFrame(0.,0.,1.1,4.0);

  hFrame->SetTitle(";p_{T}^{MC} (GeV/c); DCA resolution, cm");
  hFrame->SetTitle(Form("#eta = %.1f",eta));
  hFrame->GetXaxis()->SetTitleOffset(1.2);
  hFrame->GetXaxis()->SetTitleSize(0.045);
  hFrame->GetXaxis()->SetLabelSize(0.045);
  hFrame->GetYaxis()->SetTitleSize(0.045);
  hFrame->GetYaxis()->SetLabelSize(0.045);

  gRes1->SetLineStyle(1);
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
  
  gPad->Print(Form("dca_resolution%s%.2f.png",pid.Data(),eta));
}

