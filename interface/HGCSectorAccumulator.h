#ifndef _hgcsectoraccumulator_h_
#define _hgcsectoraccumulator_h_

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DetectorDescription/Base/interface/DDTranslation.h"
#include "DetectorDescription/Base/interface/DDRotationMatrix.h"

#include "TMath.h"
#include "TString.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TVector2.h"

/**
   @class HGCSectorAccumulator
   @short wrapper for sector histogramming
 */
class HGCSectorAccumulator
{
 public:
  
  HGCSectorAccumulator(int subDet,int layer, int copy);
  
  inline void setGeometry(float h, float bl,float tl, float cell, float recoCell)    
    { 
      h_=h;  bl_=bl; tl_=tl; cell_=cell; 
    }
  
  inline void setRotation(DD3Vector &xrot,DD3Vector &yrot,DD3Vector &zrot) 
    {
      xx_=xrot.x(); xy_=yrot.x(); xz_=zrot.x();
      yx_=xrot.y(); yy_=yrot.y(); yz_=zrot.y();
      zx_=xrot.z(); zy_=yrot.z(); zz_=zrot.z();
    }
  
  inline void setTranslation(DDTranslation &transl)
    {
      gx_=transl.x(); gy_=transl.y(); gz_=transl.z();
    }

  inline void setBasePhi(float basePhi)
    {
      basePhi_=basePhi;
    }
  inline float getBasePhi() { return basePhi_; }

  void configure(edm::Service<TFileService> &fs);
  int acquire(float edep, float t, float x, float y);         
  int digitize(float adc, float x, float y);         
  void reset();

  TVector3 getGlobalPointAt(int bin);
  TVector2 getLocalPointAt(int bin);
  float getEnergyDepAt(int bin);
  float getAverageTimeAt(int bin);
  TVector3 getRecoGlobalPointAt(int bin);
  TVector2 getRecoLocalPointAt(int bin);
  float getADCsAt(int bin);
  float getCellSize()     { return cell_;     }
  float getRecoCellSize() { return recoCell_; }
  float getHalfHeight()   { return h_;        }
  float getHalfBottom()   { return bl_;       }
  float getHalfTop()      { return tl_;       }
  TH2F *getAccumulator()  { return edepH_;    }
  TH2F *getTimer()        { return tH_;       }
  TH2F *getDigis()        { return adcH_;     }

  void dumpGeometry()
  {
    std::cout << title_ << std::endl
	      << "h=" << h_ << " b=" << bl_ << " t=" << tl_ << " cell=" << cell_ << std::endl
	      << "(x0,y0,z0)=(" << gx_ << "," << gy_ << "," << gz_ << ")" << std::endl
	      << "rho=" << sqrt(gx_*gx_+gy_*gy_) << "phi=" << basePhi_ << std::endl
	      << "|xx,xy,xz|=|" << xx_ << "\t" << xy_ << "\t" << xz_ << "|" << std::endl
	      << "|yx,yy,yz|=|" << yx_ << "\t" << yy_ << "\t" << yz_ << "|" << std::endl
	      << "|zx,zy,zz|=|" << zx_ << "\t" << zy_ << "\t" << zz_ << "|" << std::endl;
  }

  ~HGCSectorAccumulator();

 private:

  char id_[25], title_[50];
  float h_, bl_, tl_, cell_, recoCell_;
  float gx_, gy_, gz_, basePhi_;
  float xx_,xy_,xz_,yx_,yy_,yz_,zx_,zy_,zz_;

  TH2F *gxH_,    *gyH_,    *gzH_,    *edepH_, *tH_;
  TH2F *gxRecoH_,*gyRecoH_,*gzRecoH_,*adcH_;
};


#endif
