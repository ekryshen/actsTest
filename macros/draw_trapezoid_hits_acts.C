#if !defined(__CINT__) && !defined(__CLING__)
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "MpdFtdPoint.h"
#include "MpdFtdGeo.h"
#endif

void draw_trapezoid_hits_acts(){
  MpdFtdGeo ftdGeo;

  TFile* fHits = new TFile("acts/hits.root");
  TTree* tHits = (TTree*) fHits->Get("hits");
  float tz = 0;
  float tx = 0;
  float ty = 0;
  UInt_t event_id = 0;
  ULong64_t geometry_id = 0;
  tHits->SetBranchAddress("event_id",&event_id);
  tHits->SetBranchAddress("geometry_id",&geometry_id);
  tHits->SetBranchAddress("tx",&tx);
  tHits->SetBranchAddress("ty",&ty);
  tHits->SetBranchAddress("tz",&tz);

  TFile* fMeas = new TFile("acts/measurements.root");
  TTree* tMeas = (TTree*) fMeas->Get("measurements");
  int32_t m_event_id;
  int32_t m_volume_id;
  int32_t m_layer_id;
  int32_t m_surface_id;
  float m_rec_loc0;
  float m_rec_loc1;
  float m_true_loc0;
  float m_true_loc1;
  float m_true_z;
  tMeas->SetBranchAddress("event_nr",&m_event_id);
  tMeas->SetBranchAddress("volume_id",&m_volume_id);
  tMeas->SetBranchAddress("layer_id",&m_layer_id);
  tMeas->SetBranchAddress("surface_id",&m_surface_id);
  tMeas->SetBranchAddress("rec_loc0",&m_rec_loc0);
  tMeas->SetBranchAddress("rec_loc1",&m_rec_loc1);
  tMeas->SetBranchAddress("true_loc0",&m_true_loc0);
  tMeas->SetBranchAddress("true_loc1",&m_true_loc1);
  tMeas->SetBranchAddress("true_z",&m_true_z);

  TH1D* h = new TH1D("h","h",1000,-100,100);
  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    if (ftdGeo.GetLayerType(m_layer_id-2)>=7) continue;
    h->Fill(m_rec_loc0 - m_true_loc0);
  }
  new TCanvas;
  h->Draw();
}
