#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
//#include "MpdMCTrack.h"
//#include "MpdFwdPoint.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"

void draw_measurements(TString dir1 = "roc_pi_16_7deg", TString dir2 = "roc_pi_19_7deg"){
  gStyle->SetOptStat(0);
  dir1.Append("/");
  dir2.Append("/");

  TFile* f1 = new TFile(dir1 + "tracking_efficiency.root");
  TFile* f2 = new TFile(dir2 + "tracking_efficiency.root");
  TH1D* hMeasurements16 = (TH1D*) f1->Get("hMeasurements16");
  TH1D* hMeasurements19 = (TH1D*) f2->Get("hMeasurements19");

  TCanvas* c = new TCanvas("c", "c", 1200, 900);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.06);
  hMeasurements19->SetTitle(";Number of measurements per track;Entries");
  hMeasurements19->Draw();
  hMeasurements16->Draw("same");
  hMeasurements19->SetLineColor(kMagenta);
  hMeasurements16->SetLineColor(kBlue);
  hMeasurements19->SetLineWidth(2);
  hMeasurements16->SetLineWidth(2);
  TLegend* legend = new TLegend(0.75,0.68,0.95,0.85);
  legend->SetBorderSize(0);
  legend->AddEntry(hMeasurements19,"#eta = 1.9","l");
  legend->AddEntry(hMeasurements16,"#eta = 1.6","l");
  legend->Draw();
  
  gPad->Print(dir1 + "measurements.png");
}

