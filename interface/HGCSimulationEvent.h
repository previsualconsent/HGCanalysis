#ifndef _hgcsimulationevent_h_
#define _hgcsimulationevent_h_

#define MAXHGCHITSPEREVENT 1000000

typedef struct { 
  Int_t event, lumi, run;
  Int_t ngen, gen_id[MAXHGCHITSPEREVENT];
  Float_t gen_pt[MAXHGCHITSPEREVENT], gen_eta[MAXHGCHITSPEREVENT], gen_phi[MAXHGCHITSPEREVENT], gen_en[MAXHGCHITSPEREVENT];

  Int_t nhits;
  Int_t hit_type[MAXHGCHITSPEREVENT], hit_layer[MAXHGCHITSPEREVENT], hit_sec[MAXHGCHITSPEREVENT], hit_bin[MAXHGCHITSPEREVENT];
  Float_t hit_edep[MAXHGCHITSPEREVENT],hit_avgt[MAXHGCHITSPEREVENT];

  Int_t ndigis;
  Int_t digi_type[MAXHGCHITSPEREVENT], digi_layer[MAXHGCHITSPEREVENT], digi_sec[MAXHGCHITSPEREVENT], digi_bin[MAXHGCHITSPEREVENT];
  Float_t digi_adc[MAXHGCHITSPEREVENT];

} HGCSimEvent_t;



#endif
