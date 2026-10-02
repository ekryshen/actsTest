#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "style.h"

void analyse_dca(TString dir = ".", TString pid = "Pi", double etaMean = 1.60, bool refit = 0){
  dir.Append("/");
//  gStyle->SetStatFontSize(0.08);
  gStyle->SetStatH(0.15);
  gStyle->SetStatW(0.2);
  gStyle->SetStatFormat("6.3g");
  gStyle->SetOptFit(101);

//  TFile* f = new TFile(dir + (refit ? Form("tracking_performance_%.2f_refit.root",etaMean) : Form("tracking_performance_%.2f.root",etaMean)));
  TFile* f = new TFile(Form("tracking_performance_%.2f.root",etaMean));
  f->ls();
  
  TH2D* hResDvsPtMC = (TH2D*) f->Get(Form("hResDvsPt%s",pid.Data()));
  int nbins = hResDvsPtMC->GetNbinsX();
  printf("%d\n",nbins);
  TCanvas* c1 = new TCanvas("c1","c1",1800,800);
  c1->Divide(4,2,0.001,0.02);
  
  
  const int nPtBins = 18;
  double vPt[nPtBins];
  double vRes[nPtBins];

  TF1* fGaus = new TF1("fGaus","gaus", -1, 1);
  int imin = 3;
  for (int ibin = imin; ibin<=nbins; ibin++){
    double minPtMC = 1000*hResDvsPtMC->GetXaxis()->GetBinLowEdge(ibin);
    double maxPtMC = 1000*hResDvsPtMC->GetXaxis()->GetBinUpEdge(ibin);
    TH1D* hProj = hResDvsPtMC->ProjectionY(Form("hRes_%.0f_%.0f",minPtMC,maxPtMC),ibin,ibin);
    double mean = hProj->GetMean();
    double sigma = hProj->GetRMS();

    hProj->Fit(fGaus,"LQ0","",mean-2*sigma,mean+2*sigma);
    vPt[ibin-imin] = hResDvsPtMC->GetXaxis()->GetBinCenter(ibin);
    vRes[ibin-imin] = fGaus->GetParameter(2);
    //vRes[ibin-imin] = hProj->GetRMS();
    if (ibin%2==1) continue;
    c1->cd(ibin/2-1);
    SetPad(gPad);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.15);
    SetHisto(hProj,Form("p_{T}: %.0f - %.0f MeV;DCA resolution",minPtMC,maxPtMC));
    hProj->SetTitleOffset(1.3);
//    hProj->GetXaxis()->SetRangeUser(-0.7,0.7);
    hProj->Draw();
    hProj->GetListOfFunctions()->At(0)->Draw("same");
    
  }
  c1->Print(dir + Form("res%s%.0f_%d.png",pid.Data(), etaMean*10,int(refit)));

  for (int i=0;i<nPtBins;i++){
    printf("%d %f %f\n",i,vPt[i],vRes[i]);
  }
  TGraph* g = new TGraph(nPtBins,vPt,vRes);
  TFile* fg = new TFile(dir + (refit ? "dca_resolution_refit.root" : "dca_resolution.root"),"update");
//  TFile* fg = new TFile(dir + (refit ? "resolution_refit.root" : "mom_resolution.root"),"update");
  g->Write(Form("gDcaRes%s%.2f",pid.Data(), etaMean));
  fg->Close();
}
