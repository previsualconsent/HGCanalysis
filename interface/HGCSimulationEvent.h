#ifndef _hgcsimulationevent_h_
#define _hgcsimulationevent_h_

#define MAXHITS 1000000

typedef struct 
{
  double halfHeight, halfBottom, halfTop, halfWidth;
  double globalX,globalY,globalZ, basePhi;
  double xx,xy,xz,yx,yy,yz,zx,zy,zz;
} SectorGeometry_t;


typedef struct { 
  Int_t event, lumi, run;
  Int_t ngen, gen_id[MAXHITS];
  Float_t gen_pt[MAXHITS], gen_eta[MAXHITS], gen_phi[MAXHITS], gen_en[MAXHITS];
  Int_t nhits, hit_type[MAXHITS], hit_zp[MAXHITS], hit_layer[MAXHITS], hit_sec[MAXHITS], hit_subsec[MAXHITS], hit_cell[MAXHITS];
  Float_t hit_edep[MAXHITS], hit_x[MAXHITS], hit_y[MAXHITS], hit_t[MAXHITS], hit_gx[MAXHITS], hit_gy[MAXHITS], hit_gz[MAXHITS];
} HGCSimEvent_t;



#endif
