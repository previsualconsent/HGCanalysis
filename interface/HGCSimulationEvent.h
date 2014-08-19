#ifndef _hgcsimulationevent_h_
#define _hgcsimulationevent_h_

#include "TTree.h"

#define MAXGENPEREVENT 10
#define MAXG4PEREVENT 1000
#define MAXTKSPEREVENT 1000
#define MAXHGCHITSPEREVENT 1000000

typedef struct { 

  Int_t run, lumi, event;

  //generator level
  Int_t ngen, gen_id[MAXHGCHITSPEREVENT];
  Float_t gen_pt[MAXHGCHITSPEREVENT], gen_eta[MAXHGCHITSPEREVENT], gen_phi[MAXHGCHITSPEREVENT], gen_en[MAXHGCHITSPEREVENT];
  
  //geant4 information
  Int_t   ng4,g4_id[MAXG4PEREVENT];
  Float_t g4_vtx[MAXG4PEREVENT],g4_vty[MAXG4PEREVENT],g4_vtz[MAXG4PEREVENT];
  Float_t g4_en[MAXG4PEREVENT], g4_eta[MAXG4PEREVENT],g4_phi[MAXG4PEREVENT];
  Float_t dEnInTracker;
  
  //tracks and their extrapolation to HGCal
  Int_t ntk,tk_nhits[MAXTKSPEREVENT];
  Float_t tk_pt[MAXTKSPEREVENT],tk_eta[MAXTKSPEREVENT],tk_phi[MAXTKSPEREVENT],tk_chi2[MAXTKSPEREVENT];
  Int_t tk_extrapol_sd[MAXTKSPEREVENT][55],tk_extrapol_layer[MAXTKSPEREVENT][55];
  Float_t tk_extrapol_x[MAXTKSPEREVENT][55],tk_extrapol_y[MAXTKSPEREVENT][55];

  //sim hits and ADC counts
  Int_t nhits;
  Int_t hit_type[MAXHGCHITSPEREVENT], hit_layer[MAXHGCHITSPEREVENT];
  Float_t hit_x[MAXHGCHITSPEREVENT],hit_y[MAXHGCHITSPEREVENT],hit_z[MAXHGCHITSPEREVENT];
  Float_t hit_eta[MAXHGCHITSPEREVENT],hit_phi[MAXHGCHITSPEREVENT];
  Float_t hit_edep[MAXHGCHITSPEREVENT],hit_edep_sample[MAXHGCHITSPEREVENT][9];
  Short_t digi_adc[MAXHGCHITSPEREVENT];
  
} HGCSimEvent_t;


void initHGCSimulationEventTree(TTree *t,HGCSimEvent_t &simEvt);


#endif
