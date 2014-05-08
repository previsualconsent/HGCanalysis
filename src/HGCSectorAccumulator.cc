#include "UserCode/HGCanalysis/interface/HGCSectorAccumulator.h"

//
HGCSectorAccumulator::HGCSectorAccumulator(int subDet,int layer, int copy) : gxH_(0), gyH_(0), gzH_(0), edepH_(0)
{
  //identifier for the histos
  sprintf(id_,"s_%d_%d_%d",subDet,layer,copy);
  sprintf(title_,"Layer %d Sector %d",layer,copy);
}


//
HGCSectorAccumulator::~HGCSectorAccumulator()
{
}

//
void HGCSectorAccumulator::configure(edm::Service<TFileService> &fs)
{
  int ndivx=TMath::Floor(bl_/cell_);
  ndivx=(ndivx+TMath::Floor((tl_-ndivx*cell_)/cell_));
  int ndivy=TMath::Floor(h_/cell_);
  
  //init histos
  edepH_ = fs->make<TH2F>(TString("E_") +id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-cell_*ndivx,cell_*ndivx,2*ndivy,-cell_*ndivy,cell_*ndivy);
  gxH_   = fs->make<TH2F>(TString("gx_")+id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-cell_*ndivx,cell_*ndivx,2*ndivy,-cell_*ndivy,cell_*ndivy);
  gyH_   = fs->make<TH2F>(TString("gy_")+id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-cell_*ndivx,cell_*ndivx,2*ndivy,-cell_*ndivy,cell_*ndivy);
  gzH_   = fs->make<TH2F>(TString("gz_")+id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-cell_*ndivx,cell_*ndivx,2*ndivy,-cell_*ndivy,cell_*ndivy);

  //cache global quantities
  float rho=sqrt(pow(gx_,2)+pow(gy_,2));
  for(int xbin=1; xbin<gxH_->GetXaxis()->GetNbins(); xbin++)
    {
      for(int ybin=1; ybin<gxH_->GetYaxis()->GetNbins(); ybin++)
	{
	  float localX=gxH_->GetXaxis()->GetBinCenter(xbin);
	  float rotLocalX=localX;	
	  float localY=gxH_->GetYaxis()->GetBinCenter(ybin);
  	  float rotLocalY=localY+rho;

	  float gx = rotLocalX*(TMath::Cos(basePhi_)*xx_-TMath::Sin(basePhi_)*yx_) + rotLocalY*(TMath::Cos(basePhi_)*xy_-TMath::Sin(basePhi_)*yy_);
	  gxH_->SetBinContent(xbin,ybin,gx);

	  float gy = rotLocalX*(TMath::Sin(basePhi_)*xx_+TMath::Cos(basePhi_)*yx_) + rotLocalY*(TMath::Sin(basePhi_)*xy_+TMath::Cos(basePhi_)*yy_);
	  gyH_->SetBinContent(xbin,ybin,gy);

	  float gz = gz_;
	  gzH_->SetBinContent(xbin,ybin,gz);
	}
    }
}

//
int HGCSectorAccumulator::acquire(float edep, float t, float x, float y)
{
  if(edepH_==0) return -1;
  return edepH_->Fill(x,y,edep);
}

//
void HGCSectorAccumulator::reset()
{
  if(edepH_) edepH_->Reset("ICE");
}

//
TVector3 HGCSectorAccumulator::getGlobalPointAt(int bin)
{
  TVector3 xyz(gxH_->GetBinContent(bin),
	       gyH_->GetBinContent(bin),
	       gzH_->GetBinContent(bin) 
	       );
  return xyz;
}

//
TVector2 HGCSectorAccumulator::getLocalPointAt(int bin)
{
  Int_t binx,biny,binz;
  gxH_->GetBinXYZ(bin,binx,biny,binz);
  TVector2 xy(gxH_->GetXaxis()->GetBinCenter(binx), gxH_->GetYaxis()->GetBinCenter(biny));
  return xy;
}

//
float HGCSectorAccumulator::getEnergyDepAt(int bin)
{
  return edepH_->GetBinContent(bin);
}
