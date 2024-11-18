#include "TFile.h"
#include "TTree.h"
double rMinStation = 35.7*10; // mm
double rMaxStation = 130.*10; // mm


void analyse_occupancy(TString dir = "../macros_backup/mpd"){
//void analyse_occupancy(TString dir = "../macros_backup/mpd_noframe/"){
//void analyse_occupancy(TString dir = "../macros_backup/mpd_notpc/"){
  dir.Append("/");
  TFile* fHits = new TFile(dir + "hits.root");
  TTree* tHits = (TTree*) fHits->Get("hits");

  TFile* fPart = new TFile(dir + "particles.root");
  TTree* tPart = (TTree*) fPart->Get("particles");
  UInt_t part_event_id;
  auto part_id  = new std::vector<unsigned long>; 
  tPart->SetBranchAddress("event_id",&part_event_id);
  tPart->SetBranchAddress("particle_id",&part_id);
  tPart->BuildIndex("event_id");

  int nEvents = tPart->GetEntries();
  printf("nEvents = %d\n",nEvents);
  float tz = 0;
  float tx = 0;
  float ty = 0;
  UInt_t event_id = 0;
  ULong64_t particle_id = 0;
  ULong64_t geometry_id = 0;
  tHits->SetBranchAddress("event_id",&event_id);
  tHits->SetBranchAddress("particle_id",&particle_id);
  tHits->SetBranchAddress("tx",&tx);
  tHits->SetBranchAddress("ty",&ty);
  tHits->SetBranchAddress("tz",&tz);

  TH1D* hZ = new TH1D("hZ","",100,2000,3100);
  TH1D* hR = new TH1D("hR","",130,0,130);
  TH1D* hRscaled = new TH1D("hRscaled","",130,0,130);

  TH1D* hMult = new TH1D("hMult","",100,0,3000);
  vector<int> vMult(tPart->GetEntries(),0);

  int nSelected = 0;
  float minMult = 2500;
  for (int ev=0;ev<tPart->GetEntries();ev++){
    tPart->GetEntry(ev);
    float mult = part_id->size(); 
    hMult->Fill(mult);
    vMult[ev] = mult;
    if (mult<minMult) continue;
    nSelected++;
  }
  printf("%d\n",nSelected);
  
  for (int entry=0;entry< tHits->GetEntries();entry++){
    tHits->GetEntry(entry);
    printf("%d\n",event_id);
    if (vMult[event_id]<minMult) continue;
    double r = sqrt(tx*tx+ty*ty);
    if (r<rMinStation || r > rMaxStation) continue;
    hZ->Fill(tz);
    if (tz>2200) continue;
    double rcm = 0.1*r;
    hR->Fill(rcm);
    hRscaled->Fill(rcm,1./2./3.14/rcm);
  }
  hR->Scale(1./nSelected);
  hRscaled->Scale(1./nSelected);
  new TCanvas;
  hZ->Draw();
  new TCanvas;
  hR->Draw();
  new TCanvas;
  hRscaled->Draw();
}

