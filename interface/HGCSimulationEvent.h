#ifndef _hgcsimulationevent_h_
#define _hgcsimulationevent_h_

#define MAXHITS 1000000

typedef struct { 
  Int_t event, lumi, run;
  Int_t ngen, gen_id[MAXHITS];
  Float_t gen_pt[MAXHITS], gen_eta[MAXHITS], gen_phi[MAXHITS], gen_en[MAXHITS];
  Int_t nee, ee_zp[MAXHITS], ee_layer[MAXHITS], ee_sec[MAXHITS], ee_subsec[MAXHITS], ee_cell[MAXHITS];
  Float_t ee_edep[MAXHITS], ee_x[MAXHITS], ee_y[MAXHITS], ee_t[MAXHITS], ee_gx[MAXHITS], ee_gy[MAXHITS], ee_gz[MAXHITS];
} HGCSimEvent_t;



#endif
