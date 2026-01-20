#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "vector"
#include "style.h"
#include "tree_summary.C"
#include "tree_trackstates.C"

void analyse_ambiguity(TString dir = "../build/test/"){
  TFile* fStates = new TFile(dir + "trackstates.root");
  TTree* tStates = (TTree*) fStates->Get("trackstates");
  SetStatesBranchAddresses(tStates);

  TFile* fTracks = new TFile(dir + "tracksummary.root");
  TTree* tTracks = (TTree*) fTracks->Get("tracksummary");
  SetBranchAddresses(tTracks);
  
  int nEvents = tStates->GetEntries();
  for (int ev=0; ev<nEvents; ev++){
    tStates->GetEntry(ev);
    tTracks->GetEntry(ev);

    for (int is=0; is<m_layerID->size(); is++){
      if (m_particleId->at(is).size()==0 || m_particleId->at(is).size()>1) continue;
      auto id = m_particleId->at(is)[0][2];
    }
  }

}
