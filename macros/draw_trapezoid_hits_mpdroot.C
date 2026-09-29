#if !defined(__CINT__) && !defined(__CLING__)
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "MpdFtdPoint.h"
#include "MpdFtdGeo.h"
#endif
#include "draw_trapezoid_station.C"

void draw_trapezoid_hits_mpdroot(int selected_event = 0, int station = 0){
  draw_trapezoid_station(0);
  
  MpdFtdGeo fg;

  auto fpoints = new TFile("mc.root");
  auto tpoints = (TTree*) fpoints->Get("mpdsim");
  auto points = new TClonesArray("MpdFtdPoint");
  tpoints->SetBranchAddress("FtdPoint",&points);

  TFile* fhits = new TFile("ftd.root");
  auto thits = (TTree*) fhits->Get("mpdsim");
  auto hits = new TClonesArray("MpdFtdHit");
  thits->SetBranchAddress("FtdHit",&hits);

  int nEvents = tpoints->GetEntries();
//  TH1D* h = new TH1D("h","h",200,-0.1,0.1);
  TH1D* h = new TH1D("h","h",200,-0.1,0.1);

  TGraph* g = new TGraph();
  g->SetMarkerColor(kMagenta);
  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerSize(1.);

  for (int ev=0;ev<nEvents;ev++){
    if (selected_event>=0 && ev!=selected_event) continue;
    tpoints->GetEntry(ev);
    thits->GetEntry(ev);

    int nPoints = points->GetEntriesFast();
    printf("nPoints = %d\n",nPoints);
    for (int ip=0;ip<nPoints;ip++){
      MpdFtdPoint* point = (MpdFtdPoint*) points->At(ip);
      int detID = point->GetDetectorID();
      int layer = detID%100;
      int module = detID/100;
      if (fg.GetLayerType(layer)<7) continue;
      if (fg.GetLayerStation(layer)!=station) continue;
      double x = point->GetX();
      double y = point->GetY();
      double z = point->GetZ();
      printf("layer=%2d module=%2d x=%6.2f y=%6.2f z=%6.2f\n", layer, module, x, y, z);
      g->AddPoint(x,y);
    }

    int nHits = hits->GetEntriesFast();
    printf("nHits = %d\n",nHits);
    for (int ih=0;ih<nHits;ih++){
      MpdFtdHit* hit = (MpdFtdHit*) hits->At(ih);
      int detID = hit->GetDetectorID();
      double loc0 = hit->GetX();
      double loc1 = hit->GetY();
      double z = hit->GetZ();
      int ip = hit->GetRefIndex();
      MpdFtdPoint* point = (MpdFtdPoint*) points->At(ip);
      double x = point->GetX();
      double y = point->GetY();
      int layer = detID%100;
      int module = detID/100;
      if (fg.GetLayerType(layer)<7) continue;
      if (fg.GetLayerStation(layer)!=station) continue;
      auto color = kOrange+7;
      if (fg.GetLayerType(layer) == 8) color = kGreen+1;
      if (fg.GetLayerType(layer) == 9) color = kAzure;
      auto [true_loc0, true_loc1] = fg.GetModuleLocCoordinates(x, y, layer, module);
      auto [xh, yh] = fg.GetModuleGlobalCoordinates(loc0, loc1, layer, module);
      printf("layer=%2d module=%2d loc0=%6.2f loc1=%6.2f x=%6.2f y=%6.2f z=%6.2f\n", layer, module, loc0, loc1, xh, yh, z);
      // h->Fill(xh-x);
      h->Fill(loc0-true_loc0);      

      auto polyLine = fg.GetHitPolyLine(true_loc0, layer, module);
      polyLine->SetLineColor(color);
      polyLine->SetLineWidth(2);
      polyLine->Draw();
    }
  }

  g->Draw("p");
  gPad->Print("trapezoid_hits.png");

  new TCanvas;
  gPad->SetLogy();
  h->Draw();
  

}
