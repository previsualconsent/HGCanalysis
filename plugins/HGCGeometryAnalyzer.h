#ifndef _HGCGeometryAnalyzer_h_
#define _HGCGeometryAnalyzer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "SimG4CMS/Calo/interface/HGCNumberingScheme.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include "TH1F.h"
#include "TH2F.h"

#include <string>

/**
   @class HGCGeometryAnalyzer
   @author P. Silva (CERN)
*/

class HGCGeometryAnalyzer : public edm::EDAnalyzer 
{
  
 public:
  
  explicit HGCGeometryAnalyzer( const edm::ParameterSet& );
  ~HGCGeometryAnalyzer();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );

 private:

  bool testDone_;
  std::vector<std::string> geometrySource_;
  edm::Service<TFileService> *fs_;
};
 

#endif
