#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "vector"
#include "style.h"
#include "tree_trackstates.C"

void analyse_trackstates(){
  // TFile* f = new TFile("../build/trackstates_ckf.root");
  TFile* f = new TFile("../build/test/trackstates.root");
  TTree* t = (TTree*) f->Get("trackstates");
  SetBranchAddresses(t);
  int nEvents = t->GetEntries();
  TH1D* hLayers = new TH1D("hLayers","",40,0,40);
  TH1D* hMeas = new TH1D("hMeas","",40,0,40);
  TH1D* hMeasLost = new TH1D("hMeasLost","",40,0,40);
  for (int ev=0; ev<nEvents; ev++){
    t->GetEntry(ev);
    printf("event = %d\n",m_eventNr);
    printf("track = %d\n",m_trackNr);
    int nMeas = 0;
    int id0 = -1;
    for (int is=0; is<m_layerID->size(); is++){
      auto type = m_stateType->at(is);
      if (type!=0) continue;
      auto x = m_t_x->at(is);
      auto y = m_t_y->at(is);
      auto layer = m_layerID->at(is);
      //printf("%d %d %f %f\n",type, layer,x ,y);
      //printf("%d\n",)
      hLayers->Fill(layer);
      if (m_particleId->at(is).size()==0) {
        printf("particle Id is empty\n");
        continue;
      } 
      if (m_particleId->at(is).size()>1) {
        printf("particle Id array is larger than 1\n");
        continue;
      }
      auto id = m_particleId->at(is)[0][2];
      printf("%d %d\n", layer, id);
      nMeas++;
      if (id0==-1) {
        id0 = id;
        continue;
      }
      else if (id0!=id and nMeas==2) {
        hMeasLost->Fill(layer);
      }
    }
    hMeas->Fill(nMeas);
  }
  hLayers->Draw();
  new TCanvas;
  hMeas->Draw();
  new TCanvas;
  hMeasLost->Draw();
}