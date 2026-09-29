#if !defined(__CINT__) && !defined(__CLING__)
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "MpdFtdPoint.h"
#include "MpdFtdGeo.h"
#endif

void draw_trapezoid_station(int station=0, bool zoom=0, double xmin=60, double ymin=0, double xmax=100, double ymax=40){
  MpdFtdGeo fg;
  int nLayersPerStation = fg.GetNLayersPerStation();
  int layerMin = station*nLayersPerStation;
  int layerMax = layerMin + nLayersPerStation - 1;
  double rmin = fg.GetLayerRMin(layerMin);
  double rmax = fg.GetLayerRMax(layerMax);  

  new TCanvas("cmeas","cmeas",850,850);
  gPad->SetRightMargin(0.005);
  gPad->SetTopMargin(0.005);
  gPad->SetLeftMargin(0.1);
  gPad->SetBottomMargin(0.07);
  if (!zoom) {
    xmin = -rmax;
    ymin = -rmax;
    xmax = +rmax;
    ymax = +rmax;
  }
  TH1F* hFrame = gPad->DrawFrame(xmin, ymin, xmax, ymax);
  hFrame->GetYaxis()->SetTitleOffset(1.5);
  hFrame->SetTitle(";x (cm); y (cm)"); 

  TEllipse* elDot = new TEllipse(0,0, 2);
  TEllipse* elMin = new TEllipse(0,0,rmin);
  TEllipse* elMax = new TEllipse(0,0,rmax);
  elDot->SetFillStyle(0);
  elMin->SetFillStyle(1001);
  elMax->SetFillStyle(1001);
  elMax->SetFillColor(kOrange-5);
  elMin->SetFillColor(kWhite);
  elMax->Draw("f");
  elMin->Draw("f");
  elDot->Draw();

  int layer = layerMin;
  int nModules = fg.GetLayerNumberOfModules(layer);
  for (int module=0; module < nModules; module++) {
    auto polyLine =  fg.GetModulePolyLine(layer, module);
    polyLine->SetFillColor(kWhite);
    polyLine->SetLineColor(kBlack);
    polyLine->SetLineWidth(2);
    polyLine->Draw("f");
    polyLine->Draw("");    
  }
}
