#ifndef _hgcsimulationevent_h_
#define _hgcsimulationevent_h_

#define MAXHGCSIMHITS 1000000

typedef struct { 
  Int_t event, lumi, run;
  Int_t ngen, gen_id[MAXHGCSIMHITS];
  Float_t gen_pt[MAXHGCSIMHITS], gen_eta[MAXHGCSIMHITS], gen_phi[MAXHGCSIMHITS], gen_en[MAXHGCSIMHITS];
  Int_t nhits, hit_type[MAXHGCSIMHITS], hit_layer[MAXHGCSIMHITS], hit_sec[MAXHGCSIMHITS], hit_bin[MAXHGCSIMHITS];
  Float_t hit_edep[MAXHGCSIMHITS];
} HGCSimEvent_t;



#endif
