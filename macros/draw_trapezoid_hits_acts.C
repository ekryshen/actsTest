#if !defined(__CINT__) && !defined(__CLING__)
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "MpdFtdPoint.h"
#include "MpdFtdGeo.h"
#endif
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "draw_trapezoid_station.C"

void draw_trapezoid_hits_acts(int selected_event = 0, int station = 0){
  TString dir = "acts/";

  draw_trapezoid_station(0);
  MpdFtdGeo fg;

  TFile* fHits = new TFile(dir + "hits.root");
  TTree* tHits = (TTree*) fHits->Get("hits");
  float tx = 0;                tHits->SetBranchAddress("tx",&tx);
  float ty = 0;                tHits->SetBranchAddress("ty",&ty);
  float tz = 0;                tHits->SetBranchAddress("tz",&tz);
  UInt_t event_id = 0;         tHits->SetBranchAddress("event_id",&event_id);
  ULong64_t geometry_id = 0;   tHits->SetBranchAddress("geometry_id",&geometry_id);

  TFile* fMeas = new TFile(dir + "measurements.root");
  TTree* tMeas = (TTree*) fMeas->Get("measurements");
  int32_t m_event_id;         tMeas->SetBranchAddress("event_nr",&m_event_id);
  int32_t m_volume_id;        tMeas->SetBranchAddress("volume_id",&m_volume_id);
  int32_t m_layer_id;         tMeas->SetBranchAddress("layer_id",&m_layer_id);
  int32_t m_surface_id;       tMeas->SetBranchAddress("surface_id",&m_surface_id);
  float m_rec_loc0;           tMeas->SetBranchAddress("rec_loc0",&m_rec_loc0);
  float m_rec_loc1;           tMeas->SetBranchAddress("rec_loc1",&m_rec_loc1);
  float m_true_loc0;          tMeas->SetBranchAddress("true_loc0",&m_true_loc0);
  float m_true_loc1;          tMeas->SetBranchAddress("true_loc1",&m_true_loc1);
  float m_true_x;             tMeas->SetBranchAddress("true_x",&m_true_x);
  float m_true_y;             tMeas->SetBranchAddress("true_y",&m_true_y);
  float m_true_z;             tMeas->SetBranchAddress("true_z",&m_true_z);

  TFile* fSpacepoints = new TFile(dir + "spacepoints.root");
  TTree* tSpacepoints = (TTree*) fSpacepoints->Get("spacepoints");
  float sx;                   tSpacepoints->SetBranchAddress("x",&sx);
  float sy;                   tSpacepoints->SetBranchAddress("y",&sy);
  float sz;                   tSpacepoints->SetBranchAddress("z",&sz);
  UInt_t sevent_id;           tSpacepoints->SetBranchAddress("event_id",&sevent_id);
  ULong64_t sgeometry_id;     tSpacepoints->SetBranchAddress("geometry_id",&sgeometry_id);
  

  TGraph* g = new TGraph();
  TH1D* h = new TH1D("h","h",1000,-100,100);
  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    if (selected_event>=0 && m_event_id!=selected_event) continue;
    int layer = m_layer_id-2;
    int module = m_surface_id-1;
    if (fg.GetLayerStation(layer)!=station) continue;
    if (fg.GetLayerType(layer)<7) continue;
    h->Fill(m_rec_loc0 - m_true_loc0);
    g->AddPoint(m_true_x/10,m_true_y/10);
    auto color = kOrange+7;
    if (fg.GetLayerType(layer) == 8) color = kGreen+1;
    if (fg.GetLayerType(layer) == 9) color = kAzure;
    auto polyLine = fg.GetHitPolyLine(m_true_loc0/10., layer, module);
    polyLine->SetLineColor(color);
    polyLine->SetLineWidth(2);
    polyLine->Draw();
  }

  g->SetMarkerColor(kMagenta);
  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerSize(1.);
  g->Draw("p");

  TGraph* gSP = new TGraph();
  for (int is=0;is<tSpacepoints->GetEntries();is++){
    tSpacepoints->GetEntry(is);
    if (selected_event>=0 && sevent_id!=selected_event) continue;
    auto geoId = Acts::GeometryIdentifier(sgeometry_id);
    int layer = geoId.layer()-2;
    int module = geoId.sensitive()-1;
    printf("%d %d\n", layer, module);
    if (fg.GetLayerStation(layer)!=station) continue;
    if (fg.GetLayerType(layer)<7) continue;
    gSP->AddPoint(sx/10.,sy/10.);
  }
  gSP->SetMarkerColor(kBlack);
  gSP->SetMarkerStyle(kOpenCross);
  gSP->SetMarkerSize(3.);
  gSP->Draw("p");

  gPad->Print("trapezoid_hits_acts.png");

  new TCanvas;
  h->Draw();
}
