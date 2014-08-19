#include "UserCode/HGCanalysis/interface/HGCSimulationEvent.h"

void initHGCSimulationEventTree(TTree *t,HGCSimEvent_t &simEvt)
{
  t->Branch("run",       &simEvt.run,        "run/I");
  t->Branch("lumi",      &simEvt.lumi,       "lumi/I");
  t->Branch("event",     &simEvt.event,      "event/I");

  t->Branch("ngen",      &simEvt.ngen,       "ngen/I");
  t->Branch("gen_id",     simEvt.gen_id,     "gen_id[ngen]/I");
  t->Branch("gen_pt",     simEvt.gen_pt,     "gen_pt[ngen]/F");
  t->Branch("gen_eta",    simEvt.gen_eta,    "gen_eta[ngen]/F");
  t->Branch("gen_phi",    simEvt.gen_phi,    "gen_phi[ngen]/F");
  t->Branch("gen_en",     simEvt.gen_en,     "gen_en[ngen]/F"); 

  t->Branch("ng4",      &simEvt.ng4,       "ng4/I");
  t->Branch("g4_id",     simEvt.g4_id,     "g4_id[ngen]/I");
  t->Branch("g4_vtx",    simEvt.g4_vtx,    "g4_vtx[ngen]/F");
  t->Branch("g4_vty",    simEvt.g4_vty,    "g4_vty[ngen]/F");
  t->Branch("g4_vtz",    simEvt.g4_vtz,    "g4_vtz[ngen]/F");
  t->Branch("g4_en",     simEvt.g4_en,     "g4_en[ngen]/F");
  t->Branch("g4_eta",    simEvt.g4_eta,    "g4_eta[ngen]/F");
  t->Branch("g4_phi",    simEvt.g4_phi,    "g4_phi[ngen]/F");
  t->Branch("g4_dEnInTracker",    &simEvt.g4_dEnInTracker,    "g4_dEnInTracker/F");

  t->Branch("ntk",              &simEvt.ntk,               "ntk/I");
  t->Branch("tk_nhits",          simEvt.tk_nhits,          "tk_nhits[ntk]/I");
  t->Branch("tk_pt",             simEvt.tk_pt,             "tk_pt[ntk]/F");
  t->Branch("tk_eta",            simEvt.tk_eta,            "tk_eta[ntk]/F");
  t->Branch("tk_phi",            simEvt.tk_phi,            "tk_phi[ntk]/F");
  t->Branch("tk_chi2",           simEvt.tk_chi2,           "tk_chi2[ntk]/F");
  t->Branch("tk_extrapol_sd",    simEvt.tk_extrapol_sd,    "tk_extrapol_sd[ntk][55]/I");
  t->Branch("tk_extrapol_layer", simEvt.tk_extrapol_layer, "tk_extrapol_layer[ntk][55]/I");
  t->Branch("tk_extrapol_x",     simEvt.tk_extrapol_x,     "tk_extrapol_x[ntk][55]/F");
  t->Branch("tk_extrapol_y",     simEvt.tk_extrapol_y,     "tk_extrapol_y[ntk][55]/F");

  t->Branch("nhits",     &simEvt.nhits,      "nhits/I");
  t->Branch("hit_type",   simEvt.hit_type,   "hit_type[nhits]/I");
  t->Branch("hit_layer",  simEvt.hit_layer,  "hit_layer[nhits]/I");
  t->Branch("hit_edep",   simEvt.hit_edep,   "hit_edep[nhits]/F");
  t->Branch("hit_edep_sample",   simEvt.hit_edep_sample,   "hit_edep_sample[nhits][9]/F");
  t->Branch("hit_x",      simEvt.hit_x,      "hit_x[nhits]/F");
  t->Branch("hit_y",      simEvt.hit_y,      "hit_y[nhits]/F");
  t->Branch("hit_z",      simEvt.hit_z,      "hit_z[nhits]/F");
  t->Branch("hit_eta",    simEvt.hit_eta,    "hit_eta[nhits]/F");
  t->Branch("hit_phi",    simEvt.hit_phi,    "hit_phi[nhits]/F");
  t->Branch("digi_adc",   simEvt.digi_adc,   "digi_adc[nhits]/s");
}
