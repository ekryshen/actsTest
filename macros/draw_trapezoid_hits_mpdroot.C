#if !defined(__CINT__) && !defined(__CLING__)
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "MpdFtdPoint.h"
#include "MpdFtdGeo.h"
#endif

void draw_trapezoid_hits_mpdroot(){
  MpdFtdGeo ftdGeo;

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
  for (int ev=0;ev<nEvents;ev++){
    tpoints->GetEntry(ev);
    thits->GetEntry(ev);

    int nPoints = points->GetEntriesFast();
    printf("nPoints = %d\n",nPoints);
    for (int ip=0;ip<nPoints;ip++){
      MpdFtdPoint* point = (MpdFtdPoint*) points->At(ip);
      int detID = point->GetDetectorID();
      double x = point->GetX();
      double y = point->GetY();
      double z = point->GetZ();
      printf("layer=%2d module=%2d x=%6.2f y=%6.2f z=%6.2f\n", detID%100, detID/100, x, y, z);
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
      auto [true_loc0, true_loc1] = ftdGeo.GetModuleLocCoordinates(x, y, layer, module);
      auto [xh, yh] = ftdGeo.GetModuleGlobalCoordinates(loc0, loc1, layer, module);
      printf("layer=%2d module=%2d loc0=%6.2f loc1=%6.2f x=%6.2f y=%6.2f z=%6.2f\n", layer, module, loc0, loc1, xh, yh, z);
      if (ftdGeo.GetLayerType(layer)<7) continue;
      // h->Fill(xh-x);
      h->Fill(loc0-true_loc0);      
    }
  }
  new TCanvas;
  gPad->SetLogy();
  h->Draw();
}
