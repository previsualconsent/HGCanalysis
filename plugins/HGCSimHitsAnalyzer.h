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

#include "DataFormats/Candidate/interface/Candidate.h"

#include "UserCode/HGCanalysis/interface/HGCSimulationEvent.h"

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
  
  void analyzeEEHits(edm::Handle<edm::PCaloHitContainer> &caloHits, edm::Handle<edm::View<reco::Candidate> > &gen);

  bool geometryDefined_;
  std::string ddViewName_;
  std::string eeHits_, heHits_, genSource_;

  TTree *t_;
  HGCSimEvent_t simEvt_;
};
 

#endif
