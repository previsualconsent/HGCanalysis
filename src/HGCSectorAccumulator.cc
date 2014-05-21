#include "UserCode/HGCanalysis/interface/HGCSectorAccumulator.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Plane3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Geometry/Transform3D.h"

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
  //build the local -> global transformation
  const CLHEP::HepRep3x3 rotation3x3 ( xx_, yx_, zx_, xy_, yy_, zy_, xz_, yz_, zz_ );
  const CLHEP::HepRotation rotationMatrix( rotation3x3 );
  const CLHEP::Hep3Vector translationVec( gx_, gy_, gz_);
  const HepGeom::Transform3D local2globalTr(rotationMatrix,translationVec);
  
  //cache global quantities
  float rho=sqrt(pow(gx_,2)+pow(gy_,2));

  //init sim histos
  int ndivx=TMath::Floor(bl_/cell_);
  ndivx=(ndivx+TMath::Floor((tl_-ndivx*cell_)/cell_));
  int ndivy=TMath::Floor(h_/cell_);
  edepH_ = fs->make<TH2F>(TString("E_")    +id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-cell_*ndivx,cell_*ndivx,2*ndivy,-cell_*ndivy,cell_*ndivy);
  tH_    = fs->make<TH2F>(TString("AvgT_") +id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-cell_*ndivx,cell_*ndivx,2*ndivy,-cell_*ndivy,cell_*ndivy);
  gxH_   = fs->make<TH2F>(TString("gx_")   +id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-cell_*ndivx,cell_*ndivx,2*ndivy,-cell_*ndivy,cell_*ndivy);
  gyH_   = fs->make<TH2F>(TString("gy_")   +id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-cell_*ndivx,cell_*ndivx,2*ndivy,-cell_*ndivy,cell_*ndivy);
  gzH_   = fs->make<TH2F>(TString("gz_")   +id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-cell_*ndivx,cell_*ndivx,2*ndivy,-cell_*ndivy,cell_*ndivy);
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

	  const HepGeom::Point3D<float> lcoord(localX,localY,0);
	  const HepGeom::Point3D<float> gcoord( local2globalTr*lcoord );

	  std::cout << "Previous method (" <<  gx << "," << gy << "," << gz << std::endl
		    << "New method      (" <<  gcoord.x() << "," << gcoord.y() << "," << gcoord.z() << std::endl << std::endl;  
	}
    }

  //init reco histos
  ndivx      = TMath::Floor(bl_/recoCell_);
  ndivx      = (ndivx+TMath::Floor((tl_-ndivx*recoCell_)/recoCell_));
  ndivy      = TMath::Floor(h_/recoCell_);
  adcH_      = fs->make<TH2F>(TString("ADC_") +id_,title_+TString(";x [mm];y [mm]"),  2*ndivx,-recoCell_*ndivx,recoCell_*ndivx,2*ndivy,-recoCell_*ndivy,recoCell_*ndivy);
  gxRecoH_   = fs->make<TH2F>(TString("recogx_")+id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-recoCell_*ndivx,recoCell_*ndivx,2*ndivy,-recoCell_*ndivy,recoCell_*ndivy);
  gyRecoH_   = fs->make<TH2F>(TString("recogy_")+id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-recoCell_*ndivx,recoCell_*ndivx,2*ndivy,-recoCell_*ndivy,recoCell_*ndivy);
  gzRecoH_   = fs->make<TH2F>(TString("recogz_")+id_,title_+TString(";x [mm];y [mm]"),2*ndivx,-recoCell_*ndivx,recoCell_*ndivx,2*ndivy,-recoCell_*ndivy,recoCell_*ndivy);
  for(int xbin=1; xbin<gxRecoH_->GetXaxis()->GetNbins(); xbin++)
    {
      for(int ybin=1; ybin<gxRecoH_->GetYaxis()->GetNbins(); ybin++)
	{
	  float localX=gxRecoH_->GetXaxis()->GetBinCenter(xbin);
	  float rotLocalX=localX;	
	  float localY=gxRecoH_->GetYaxis()->GetBinCenter(ybin);
  	  float rotLocalY=localY+rho;

	  float gx = rotLocalX*(TMath::Cos(basePhi_)*xx_-TMath::Sin(basePhi_)*yx_) + rotLocalY*(TMath::Cos(basePhi_)*xy_-TMath::Sin(basePhi_)*yy_);
	  gxRecoH_->SetBinContent(xbin,ybin,gx);
	  
	  float gy = rotLocalX*(TMath::Sin(basePhi_)*xx_+TMath::Cos(basePhi_)*yx_) + rotLocalY*(TMath::Sin(basePhi_)*xy_+TMath::Cos(basePhi_)*yy_);
	  gyRecoH_->SetBinContent(xbin,ybin,gy);
	  
	  float gz = gz_;
	  gzRecoH_->SetBinContent(xbin,ybin,gz);
	}
    }
}

//
int HGCSectorAccumulator::acquire(float edep, float t, float x, float y)
{
  if(edepH_==0) return -1;
  if(x>tl_ || x < -tl_ || y>h_ || y<-h_)
    {
      dumpGeometry();
      std::cout << " Can't accumulate @ (" << x << " " << y << ")" <<  std::endl;
      return 0;
    }
  tH_->Fill(x,y,t*edep);
  return edepH_->Fill(x,y,edep);
}

//
int HGCSectorAccumulator::digitize(float adc, float x, float y)
{
  if(adcH_==0) return -1;
  return adcH_->Fill(x,y,adc);
}


//
void HGCSectorAccumulator::reset()
{
  if(edepH_)  edepH_->Reset("ICE");
  if(tH_)     tH_->Reset("ICE");
  if(adcH_)   adcH_->Reset("ICE");
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
TVector3 HGCSectorAccumulator::getRecoGlobalPointAt(int bin)
{
  TVector3 xyz(gxRecoH_->GetBinContent(bin),
	       gyRecoH_->GetBinContent(bin),
	       gzRecoH_->GetBinContent(bin) 
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
TVector2 HGCSectorAccumulator::getRecoLocalPointAt(int bin)
{
  Int_t binx,biny,binz;
  gxRecoH_->GetBinXYZ(bin,binx,biny,binz);
  TVector2 xy(gxRecoH_->GetXaxis()->GetBinCenter(binx), gxRecoH_->GetYaxis()->GetBinCenter(biny));
  return xy;
}

//
float HGCSectorAccumulator::getEnergyDepAt(int bin)
{
  return edepH_ ? edepH_->GetBinContent(bin) : 0.;
}

//
float HGCSectorAccumulator::getAverageTimeAt(int bin)
{
  return tH_ ? tH_->GetBinContent(bin) : 0.;
}

//
float HGCSectorAccumulator::getADCsAt(int bin)
{
  return adcH_ ? adcH_->GetBinContent(bin) : 0.;
}
