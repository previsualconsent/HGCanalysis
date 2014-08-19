#include "UserCode/HGCanalysis/plugins/HGCSimHitsAnalyzer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "SimG4CMS/Calo/interface/CaloHitID.h"

#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"

#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"

#include "DataFormats/GeometryVector/interface/Basic3DVector.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Plane3D.h"
#include "CLHEP/Geometry/Vector3D.h"

#include "TVector2.h"

#include <iostream>

using namespace std;

//
HGCSimHitsAnalyzer::HGCSimHitsAnalyzer( const edm::ParameterSet &iConfig ) : geometryDefined_(false)
{
  ddViewName_      = iConfig.getUntrackedParameter<std::string>("ddViewName", "");
  genSource_       = iConfig.getUntrackedParameter<std::string>("genSource",  "genParticles");
  hitCollections_  = iConfig.getUntrackedParameter< std::vector<std::string> >("hitCollections");
  digiCollections_ = iConfig.getUntrackedParameter< std::vector<std::string> >("digiCollections");
  sdTags_          = iConfig.getUntrackedParameter< std::vector<std::string> >("sdTags");
  addGlobalPos_    = iConfig.getUntrackedParameter< bool > ("addGlobalPos");

  edm::Service<TFileService> fs;
  fs_=&fs;
  t_=fs->make<TTree>("HGC","Event Summary");

  t_->Branch("run",       &simEvt_.run,        "run/I");
  t_->Branch("lumi",      &simEvt_.lumi,       "lumi/I");
  t_->Branch("event",     &simEvt_.event,      "event/I");
  t_->Branch("ngen",      &simEvt_.ngen,       "ngen/I");
  t_->Branch("gen_id",     simEvt_.gen_id,     "gen_id[ngen]/I");
  t_->Branch("gen_pt",     simEvt_.gen_pt,     "gen_pt[ngen]/F");
  t_->Branch("gen_eta",    simEvt_.gen_eta,    "gen_eta[ngen]/F");
  t_->Branch("gen_phi",    simEvt_.gen_phi,    "gen_phi[ngen]/F");
  t_->Branch("gen_en",     simEvt_.gen_en,     "gen_en[ngen]/F"); 
  t_->Branch("nhits",     &simEvt_.nhits,      "nhits/I");
  t_->Branch("hit_type",   simEvt_.hit_type,   "hit_type[nhits]/I");
  t_->Branch("hit_layer",  simEvt_.hit_layer,  "hit_layer[nhits]/I");
  t_->Branch("hit_edep",   simEvt_.hit_edep,   "hit_edep[nhits]/F");
  t_->Branch("hit_avgt",   simEvt_.hit_avgt,   "hit_avgt[nhits]/F");
  t_->Branch("hit_x",      simEvt_.hit_x,      "hit_x[nhits]/F");
  t_->Branch("hit_y",      simEvt_.hit_y,      "hit_y[nhits]/F");
  t_->Branch("hit_z",      simEvt_.hit_z,      "hit_z[nhits]/F");
  t_->Branch("hit_eta",    simEvt_.hit_eta,    "hit_eta[nhits]/F");
  t_->Branch("hit_phi",    simEvt_.hit_phi,    "hit_phi[nhits]/F");
  t_->Branch("digi_adc",   simEvt_.digi_adc,   "digi_adc[nhits]/s");
}

//
HGCSimHitsAnalyzer::~HGCSimHitsAnalyzer()
{
}

//
void HGCSimHitsAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  cout << "[HGCSimHitsAnalyzer::analyze] start" << endl;

  //check if geometry is defined
  if(!geometryDefined_){
    //get handle to the DD
    edm::ESTransientHandle<DDCompactView> ddViewH;
    iSetup.get<IdealGeometryRecord>().get( ddViewName_, ddViewH );
    
    //safety check
    if (!ddViewH.isValid() || !ddViewH.description()) {
      cout << "[HGCSimHitsAnalyzer::analyze] Handle for DD is not valid or does not contain any valid description" << endl;
      return;
    }    

    geometryDefined_=defineGeometry(ddViewH);
    if(geometryDefined_) { std::cout << "[HGCSimHitsAnalyzer::analyze] geometry has been defined" << endl;                      }
    else                 { std::cout << "[HGCSimHitsAnalyzer::analyze] could not define geometry from DD view" << endl; return; }
  }
  
  //event header
  simEvt_.run    = iEvent.id().run();
  simEvt_.lumi   = iEvent.luminosityBlock();
  simEvt_.event  = t_->GetEntriesFast()+1; //iEvent.id().event();

  //get gen level information
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  iEvent.getByLabel(edm::InputTag(genSource_), genParticles);
  simEvt_.ngen=0;
  for(size_t i = 0; i < genParticles->size(); ++ i)
    {
      const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[i] );
      if(p.status()!=1) continue;
      if(abs(p.pdgId())!=11 && abs(p.pdgId())!=13 && abs(p.pdgId())!=211 && abs(p.pdgId())!=130) continue;
      simEvt_.gen_id[simEvt_.ngen]=p.pdgId();
      simEvt_.gen_pt[simEvt_.ngen]=p.pt();
      simEvt_.gen_eta[simEvt_.ngen]=p.eta();
      simEvt_.gen_phi[simEvt_.ngen]=p.phi();
      simEvt_.gen_en[simEvt_.ngen]=p.energy();
      simEvt_.ngen++;
    }
  
  //hits, digis
  simEvt_.nhits=0;
  for(size_t i=0; i<hitCollections_.size(); i++)
    {
      edm::Handle<edm::PCaloHitContainer> caloHits;
      iEvent.getByLabel(edm::InputTag("g4SimHits",hitCollections_[i]),caloHits); 
      analyzeHits(i,caloHits,genParticles);

      if(digiCollections_[i].find("HE") != std::string::npos)
	{
	  edm::Handle<HGCHEDigiCollection> heDigis;
	  iEvent.getByLabel(edm::InputTag("mix",digiCollections_[i]),heDigis);
	  analyzeHEDigis(i,heDigis);
	}
      else
	{
	  edm::Handle<HGCEEDigiCollection> eeDigis;
	  iEvent.getByLabel(edm::InputTag("mix",digiCollections_[i]),eeDigis);
	  analyzeEEDigis(i,eeDigis);
	}
    }

  //dump accumulators to tree and reset
  for(size_t isd=0; isd<simHitData_.size(); isd++)
    {
      //geometry information for this sub-detector
      const HGCalDDDConstants *dddConst=numberingSchemes_[isd]->getDDDConstants();

      for(HGCSimHitDataAccumulator::iterator dit=simHitData_[isd].begin(); dit!=simHitData_[isd].end(); dit++)
	{
	  uint32_t rawDetId = dit->first;
	  uint32_t layer    = (rawDetId>>19) & 0x1F;
	  uint32_t cell     =  rawDetId      & 0xFFF;
	  uint32_t subsec   = (rawDetId>>18) & 0x1;

	  Float_t edep=dit->second[0];
	  Float_t avgt=dit->second[1];
	  if(edep>0) avgt/=edep;
	  UShort_t adc=(UShort_t)dit->second[2];
	  
	  if(edep<=0) continue;
	  if(simEvt_.nhits>=MAXHGCHITSPEREVENT) break;
	  simEvt_.hit_type [simEvt_.nhits]=isd;
	  simEvt_.hit_layer[simEvt_.nhits]=layer;
	  simEvt_.hit_edep [simEvt_.nhits]=edep;
	  simEvt_.hit_avgt [simEvt_.nhits]=avgt;
	  simEvt_.digi_adc [simEvt_.nhits]=adc;

	  //get center and rotation for this layer
	  std::vector<HGCalDDDConstants::hgtrform>::const_iterator hgtrformIt=dddConst->getFirstTrForm(true);
	  std::advance(hgtrformIt,layer-1);
	  CLHEP::Hep3Vector  h3v=hgtrformIt->h3v;
	  CLHEP::HepRotation hr=hgtrformIt->hr;


	  std::pair<float,float> localxy=dddConst->locateCell(cell,layer,subsec,true);
	  simEvt_.hit_x[simEvt_.nhits]=localxy.first;
	  simEvt_.hit_y[simEvt_.nhits]=localxy.second;
	  simEvt_.hit_z[simEvt_.nhits]=0;
	  //float rho=0;//TMath::Sqrt(TMath::power(x*x+y*y+z*z);
	  float phi=0;//TMath::ATan2(y,x);
	  float eta=0;
	  //if (rho>z) eta=0.5*TMath::Log( (rho+z)/(rho-z) );
	  //else std::cout << "HGCSimHitsAnalyser hit @ " << it->first.second << " sector#" << isec << " has rho=" <<rho << " and z=" << z << "?" << endl;
	  simEvt_.hit_eta[simEvt_.nhits]=eta;
	  simEvt_.hit_phi[simEvt_.nhits]=phi;
	
	  simEvt_.nhits++;
	}
    }
 
  //fill tree
  if(simEvt_.nhits>0)  t_->Fill();
}

//
bool HGCSimHitsAnalyzer::defineGeometry(edm::ESTransientHandle<DDCompactView> &ddViewH)
{
  if(!ddViewH.isValid()) {
    std::cout << "[HGCSimHitsAnalyzer][defineGeometry] invalid DD handle" << std::endl;
    return false;
  }

  //instantiate the relevant numbering schemes  
  const DDCompactView &cview=*ddViewH;
  HGCSimHitDataAccumulator accumulatorBase;
  for(size_t isd=0; isd<sdTags_.size(); isd++)
    {
      numberingSchemes_.push_back( new HGCNumberingScheme(cview,sdTags_[isd]) );
      simHitData_.push_back( accumulatorBase );
    }
    
  
  //all done here
  return true;
}

//
void HGCSimHitsAnalyzer::analyzeHits(size_t isd,edm::Handle<edm::PCaloHitContainer> &caloHits, edm::Handle<edm::View<reco::Candidate> > &gen)
{
  //check inputs
  if(!caloHits.isValid()) return;
  
  //analyze hits
  for(edm::PCaloHitContainer::const_iterator hit_it = caloHits->begin(); hit_it != caloHits->end(); ++hit_it) 
    {
      HGCalDetId simId(hit_it->id());
      uint32_t mySubDet(ForwardSubdetector::HGCEE);
      if(isd==1) mySubDet=ForwardSubdetector::HGCHEF;
      if(isd==2) mySubDet=ForwardSubdetector::HGCHEB;
      uint32_t id = (isd==0) ?
        (uint32_t)HGCEEDetId(ForwardSubdetector(mySubDet),simId.zside(),simId.layer(),simId.sector(),simId.subsector(),simId.cell()):
        (uint32_t)HGCHEDetId(ForwardSubdetector(mySubDet),simId.zside(),simId.layer(),simId.sector(),simId.subsector(),simId.cell());
      
      HGCSimHitDataAccumulator::iterator simHitIt=simHitData_[isd].find(id);
      if(simHitIt==simHitData_[isd].end())
        {
          HGCSimHitData baseData(3,0);
          simHitData_[isd][id]=baseData;
          simHitIt=simHitData_[isd].find(id);
        }
      (simHitIt->second)[0] += hit_it->energy();
      (simHitIt->second)[1] += hit_it->energy()*hit_it->time();
    }
}


//
void HGCSimHitsAnalyzer::analyzeHEDigis(size_t isd,edm::Handle<HGCHEDigiCollection> &heDigis)
{
  //check inputs
  if(!heDigis.isValid()) return;
  
  //analyze hits
  for(HGCHEDigiCollection::const_iterator hit_it = heDigis->begin(); hit_it != heDigis->end(); ++hit_it) 
    {
      if(hit_it->size()==0) continue;
      float rawDigi=hit_it->sample(0).raw();      
      HGCHEDetId detId(hit_it->id());
      HGCSimHitDataAccumulator::iterator simHitIt=simHitData_[isd].find(detId);
      if(simHitIt==simHitData_[isd].end())
        {
          HGCSimHitData baseData(3,0);
          simHitData_[isd][detId]=baseData;
          simHitIt=simHitData_[isd].find(detId);
        }
      (simHitIt->second)[2] = rawDigi;
    }
}

//
void HGCSimHitsAnalyzer::analyzeEEDigis(size_t isd,edm::Handle<HGCEEDigiCollection> &heDigis)
{
  //check inputs
  if(!heDigis.isValid()) return;
  
  //analyze hits
  for(HGCEEDigiCollection::const_iterator hit_it = heDigis->begin(); hit_it != heDigis->end(); ++hit_it) 
    {
      if(hit_it->size()==0) continue;
      float rawDigi=hit_it->sample(0).raw();      
      HGCEEDetId detId(hit_it->id());
      HGCSimHitDataAccumulator::iterator simHitIt=simHitData_[isd].find(detId);
      if(simHitIt==simHitData_[isd].end())
        {
          HGCSimHitData baseData(3,0);
          simHitData_[isd][detId]=baseData;
          simHitIt=simHitData_[isd].find(detId);
        }
      (simHitIt->second)[2] = rawDigi;
    }
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSimHitsAnalyzer);
