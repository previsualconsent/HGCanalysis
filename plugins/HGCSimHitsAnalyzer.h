#ifndef _HGCSimHitsAnalyzer_h_
#define _HGCSimHitsAnalyzer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimG4CMS/Calo/interface/HGCNumberingScheme.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"

#include "UserCode/HGCanalysis/interface/HGCSectorAccumulator.h"
#include "UserCode/HGCanalysis/interface/HGCSimulationEvent.h"

#include "TH1F.h"
#include "TTree.h"

#include <string>

/**
   @class HGCSimHitsAnalyzer
   @author P. Silva (CERN)
*/

class HGCSimHitsAnalyzer : public edm::EDAnalyzer 
{
  
 public:
  
  explicit HGCSimHitsAnalyzer( const edm::ParameterSet& );
  ~HGCSimHitsAnalyzer();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );

 private:

  bool defineGeometry(edm::ESTransientHandle<DDCompactView> &ddViewH);
  
  void analyzeHits(size_t isd, edm::Handle<edm::PCaloHitContainer> &caloHits, edm::Handle<edm::View<reco::Candidate> > &gen);
  void analyzeHEDigis(size_t isd,edm::Handle<HGCHEDigiCollection> &heDigis);
  void analyzeEEDigis(size_t isd,edm::Handle<HGCEEDigiCollection> &eeDigis);


  std::string ddViewName_;
  bool geometryDefined_;
  std::vector<std::string> hitCollections_, digiCollections_;
  std::vector<std::string> sdTags_;
  std::string genSource_;
  bool addGlobalPos_;

  std::vector<HGCNumberingScheme *> numberingSchemes_;

  TTree *t_;
  HGCSimEvent_t simEvt_;
  std::map< std::pair<int,int>, std::vector<HGCSectorAccumulator> > allSectors_;

  edm::Service<TFileService> *fs_;
};
 

#endif
