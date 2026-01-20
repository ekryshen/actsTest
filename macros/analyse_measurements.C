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

void analyse_measurements(){
  gStyle->SetOptStat(0);
  TFile* fMeas = new TFile("../build/test/measurements.root");
  TTree* tMeas = (TTree*) fMeas->Get("measurements");
  // tMeas->Print();
  int32_t m_event_id;
  int32_t m_volume_id;
  int32_t m_layer_id;
  float m_true_loc0;
  tMeas->SetBranchAddress("event_nr",&m_event_id);
  tMeas->SetBranchAddress("volume_id",&m_volume_id);
  tMeas->SetBranchAddress("layer_id",&m_layer_id);
  tMeas->SetBranchAddress("true_loc0",&m_true_loc0);

  TH1D* hLoc0 = new TH1D("hLoc0","",1000,-10,10);
  for (int im=0; im<tMeas->GetEntries(); im++){
    tMeas->GetEntry(im);
    if ((m_layer_id-2)%7==3) continue;
    // printf("%f\n",m_true_loc0);
    hLoc0->Fill(m_true_loc0);
  }
  new TCanvas;
  hLoc0->Draw();
}

