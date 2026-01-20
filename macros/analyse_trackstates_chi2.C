#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "vector"
#include "style.h"
#include "tree_trackstates.C"

void analyse_trackstates_chi2(){
  printf("%f\n",TMath::ChisquareQuantile(0.974,1));

  return;
  TFile* f = new TFile("../build/test/trackstates.root");
  TTree* t = (TTree*) f->Get("trackstates");
  SetBranchAddresses(t);
  int nEvents = t->GetEntries();
  TH1D* hLayers = new TH1D("hLayers","",40,0,40);
  TH1D* hMeas = new TH1D("hMeas","",40,0,40);
  TH1D* hMeasLost = new TH1D("hMeasLost","",40,0,40);
  TH2D* hLayersChi2  = new TH2D("hLayersChi2","",40,0,40,100,0,5);
  TH2D* hRes = new TH2D("hRes","",40,0,40,200,-0.5,0.5);
  TH2D* hErr = new TH2D("hErr","",40,0,40,200,-0.5,0.5);
  for (int ev=0; ev<nEvents; ev++){
    t->GetEntry(ev);
    // printf("event = %d\n",m_eventNr);
    // printf("track = %d\n",m_trackNr);
    int nMeas = 0;
    int id0 = -1;
    vector<int> particleIds;
    for (int is=0; is<m_layerID->size(); is++){
      auto type = m_stateType->at(is);
      if (type!=0) continue;
      auto x = m_t_x->at(is);
      auto y = m_t_y->at(is);
      auto layer = m_layerID->at(is);
      auto chi2 =  m_chi2->at(is);
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
      particleIds.push_back(id);
      //printf("%d %d\n", layer, id);
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
    
    if (particleIds.size()<10) continue;
    if (particleIds.size()==0) continue;
    bool isGoodTrack = 1;
    int previousId = particleIds[0];
    for (auto particleId : particleIds) {
      if (particleId!=previousId){
        isGoodTrack = 0;
        break;
      }
    }
    if (isGoodTrack == 0) continue;
    printf("good track %d: particleId=%d\n", m_trackNr, previousId);
    for (int is=0; is<m_layerID->size(); is++){
      if (m_stateType->at(is)!=0) continue; // measurements only
      auto res = m_res_eLOC0[ePredicted]->at(is);
      auto err = m_err_eLOC0[ePredicted]->at(is);
      auto layer = m_layerID->at(is);
      auto chi2 =  m_chi2->at(is);
      hLayersChi2->Fill(layer, chi2);
      hRes->Fill(layer, res);
      hErr->Fill(layer, err);
    }
  }
  hLayers->Draw();
  new TCanvas;
  hMeas->Draw();
  new TCanvas;
  hMeasLost->Draw();
  new TCanvas;
  hLayersChi2->Draw("colz");
  new TCanvas;
  int l = 4;
  TH1D* h3 = hLayersChi2->ProjectionY(Form("layer%d",l),l,l);
  h3->Draw();
  h3->Scale(1./h3->Integral());
  new TCanvas;
  auto h3cum = h3->GetCumulative();
  h3cum->Draw();
  new TCanvas;
  hRes->Draw();
  new TCanvas;
  TH1D* hResProj = hRes->ProjectionY(Form("res%d",l),l,l);
  hResProj->Draw();
  hResProj->Fit("gaus");
  new TCanvas;
  hErr->Draw();
}