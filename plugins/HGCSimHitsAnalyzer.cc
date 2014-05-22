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
  t_->Branch("hit_sec",    simEvt_.hit_sec,    "hit_sec[nhits]/I");
  t_->Branch("hit_bin",    simEvt_.hit_bin,    "hit_bin[nhits]/I");
  t_->Branch("hit_edep",   simEvt_.hit_edep,   "hit_edep[nhits]/F");
  t_->Branch("hit_avgt",   simEvt_.hit_avgt,   "hit_avgt[nhits]/F");

  t_->Branch("ndigis",    &simEvt_.ndigis,     "ndigis/I");
  t_->Branch("digi_type",  simEvt_.digi_type,  "digi_type[ndigis]/I");
  t_->Branch("digi_layer", simEvt_.digi_layer, "digi_layer[ndigis]/I");
  t_->Branch("digi_sec",   simEvt_.digi_sec,   "digi_sec[ndigis]/I");
  t_->Branch("digi_bin",   simEvt_.digi_bin,   "digi_bin[ndigis]/I");
  t_->Branch("digi_adc",   simEvt_.digi_adc,   "digi_adc[ndigis]/F");
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
      if(abs(p.pdgId())!=11 && abs(p.pdgId())!=13 && abs(p.pdgId())!=211) continue;
      simEvt_.gen_id[simEvt_.ngen]=p.pdgId();
      simEvt_.gen_pt[simEvt_.ngen]=p.pt();
      simEvt_.gen_eta[simEvt_.ngen]=p.eta();
      simEvt_.gen_phi[simEvt_.ngen]=p.phi();
      simEvt_.gen_en[simEvt_.ngen]=p.energy();
      simEvt_.ngen++;
    }
  
  //hits, digis
  simEvt_.nhits=0;
  simEvt_.ndigis=0;
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
  for(std::map< std::pair<int,int>, std::vector<HGCSectorAccumulator> >::iterator it=allSectors_.begin(); it!=allSectors_.end(); it++)
    {
      for(size_t isec=0; isec<it->second.size(); isec++)
	{
	  TH2F *accH=(it->second)[isec].getAccumulator();
	  if(accH==0) continue;
	  TH2F *tH  =(it->second)[isec].getTimer();
	  for(int xbin=1; xbin<=accH->GetXaxis()->GetNbins(); xbin++)
	    {
	      for(int ybin=1; ybin<=accH->GetYaxis()->GetNbins(); ybin++)
		{
		  float edep=accH->GetBinContent(xbin,ybin);
		  if(edep<=0) continue;
		  if(simEvt_.nhits>=MAXHGCHITSPEREVENT) break;
		  simEvt_.hit_type [simEvt_.nhits]=it->first.first;
		  simEvt_.hit_layer[simEvt_.nhits]=it->first.second;
		  simEvt_.hit_sec  [simEvt_.nhits]=isec;
		  simEvt_.hit_bin  [simEvt_.nhits]=accH->GetBin(xbin,ybin);
		  simEvt_.hit_edep [simEvt_.nhits]=edep;
		  simEvt_.hit_avgt [simEvt_.nhits]=edep > 0. ? tH->GetBinContent(xbin,ybin)/edep : 0.;
		  simEvt_.nhits++;
		}
	    }  	    
	  accH->Reset("ICE");
	  tH->Reset("ICE");

	  /*
	  TH2F *adcH=(it->second)[isec].getDigis();
	  for(int xbin=1; xbin<=adcH->GetXaxis()->GetNbins(); xbin++)
	    {
	      for(int ybin=1; ybin<=adcH->GetYaxis()->GetNbins(); ybin++)
		{
		  float adc=adcH->GetBinContent(xbin,ybin);
		  if(adc<=0) continue;
		  if(simEvt_.ndigis>=MAXHGCHITSPEREVENT) break;
		  simEvt_.digi_type [simEvt_.ndigis]=it->first.first;
		  simEvt_.digi_layer[simEvt_.ndigis]=it->first.second;
		  simEvt_.digi_sec  [simEvt_.ndigis]=isec;
		  simEvt_.digi_bin  [simEvt_.ndigis]=adcH->GetBin(xbin,ybin);
		  simEvt_.digi_adc  [simEvt_.ndigis]=adc;
		  simEvt_.ndigis++;
		}
	    }  	    
	  adcH->Reset("ICE");
	  */
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
  for(size_t isd=0; isd<sdTags_.size(); isd++)
    {
      numberingSchemes_.push_back( new HGCNumberingScheme(cview,sdTags_[isd]) );

      //instantiate energy accumulators for the different layers  
      for(std::vector<HGCalDDDConstants::hgtrform>::const_iterator trformIt=numberingSchemes_[isd]->getDDDConstants()->getFirstTrForm(); 
	  trformIt!=numberingSchemes_[isd]->getDDDConstants()->getLastTrForm(); 
	  trformIt++) 
	{
	  int layerKey(trformIt->lay);
	  std::vector<HGCalDDDConstants::hgtrap>::const_iterator modIt    = numberingSchemes_[isd]->getDDDConstants()->getFirstModule();
	  std::vector<HGCalDDDConstants::hgtrap>::const_iterator recModIt = numberingSchemes_[isd]->getDDDConstants()->getFirstModule(true);
	  for(int ilayer=1; ilayer<abs(layerKey); ilayer++) { modIt++; recModIt++; }
	  if(layerKey!=modIt->lay || layerKey!=recModIt->lay)
	    {
	      std::cout << "Inconsistency found in layer " << layerKey << " between transformation struct and geometry struct for isd= " << isd << endl;
	      continue;
	    }

	  //assign -1 if on the negative z
	  layerKey *= trformIt->zp;

	  //initiate new layer if needed
	  std::pair<int,int> sectorKey(isd,layerKey);
	  if(allSectors_.find(sectorKey)==allSectors_.end())
	    {
	      std::vector<HGCSectorAccumulator> layerSectors;
	      allSectors_[sectorKey]=layerSectors;
	    }

	  //initiate up to a new sector if needed
	  int sector=trformIt->sec;
	  if((int)allSectors_[sectorKey].size()<sector)
	    {
	      for(int isec=(int)allSectors_[sectorKey].size(); isec<sector; isec++)
		allSectors_[sectorKey].push_back( HGCSectorAccumulator(isd,layerKey,isec) );
	    }

	  //configure this specific sector
	  int isec=sector-1;
	  double simCellSize(modIt->cellSize);         //Geant4 is in mm 
	  double recoCellSize(recModIt->cellSize*10); //it comes in cm (CMS default)
	  allSectors_[sectorKey][isec].setGeometry(modIt->h, modIt->bl, modIt->tl, simCellSize,recoCellSize);

	  const HepGeom::Transform3D local2globalTr( trformIt->hr, trformIt->h3v );
	  allSectors_[sectorKey][isec].setLocal2GlobalTransformation(local2globalTr);
	  allSectors_[sectorKey][isec].configure(*fs_);    
	}
    
      for(std::map< std::pair<int,int>, std::vector<HGCSectorAccumulator> >::iterator it = allSectors_.begin();
	  it!=allSectors_.end();
	  it++)
	{
	  std::cout << it->second.size() << "sectors added for sd=" << it->first.first << " @ layer=" << it->first.second << std::endl;
	  // " with sim cell size=" << simCellSize << " and reco cell size= " << recoCellSize <<std::endl;
	}
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
      HGCalDetId detId(hit_it->id());

      int layer=detId.layer();
      int zpos=detId.zside();
      int layerKey(layer);
      if(zpos<0) layerKey *= -1;
      std::pair<int,int> sectorKey(isd,layerKey);
      int sector=detId.sector();
      //int subSector=detId.subsector();
      
      if(allSectors_.find(sectorKey)==allSectors_.end()){
	std::cout << "[HGCSimHitsAnalyzer][analyzeHits] unable to find layer parameters for detId=0x" << hex << uint32_t(detId) << dec << " iSD=" << isd << " layer=" << layerKey << std::endl;
	continue;
      }
      else if(allSectors_[sectorKey].size()<size_t(sector))
	{
	  std::cout << "[HGCSimHitsAnalyzer][analyzeHits] unable to find sector parameters for detId=0x" << hex << uint32_t(detId) << dec << " iSD=" << isd << " layer=" << layerKey << " sector=" << sector << std::endl;
	  continue;
	}

      //get local coordinates
      std::pair<float,float> xy = numberingSchemes_[isd]->getDDDConstants()->locateCell(detId.cell(),detId.layer(),detId.subsector(),false);
      float localX(xy.first), localY(xy.second);

      //accumulate energy
      //      int fbin=
      allSectors_[sectorKey][sector-1].acquire(hit_it->energy(),hit_it->time(),localX,localY);
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

      int layer=detId.layer();
      int zpos=detId.zside();
      int layerKey(layer);
      if(zpos<0) layerKey *= -1;
      std::pair<int,int> sectorKey(isd,layerKey);
      int sector=detId.sector();
      int subsector=detId.subsector();
      int cell=detId.cell();
          
      if(allSectors_.find(sectorKey)==allSectors_.end()){
	std::cout << "[HGCSimHitsAnalyzer][analyzeDigis] unable to find layer parameters for detId=0x" << hex << uint32_t(detId) << dec << " iSD=" << isd << " layer=" << layerKey << std::endl;
	continue;
      }
      else if(allSectors_[sectorKey].size()<size_t(sector))
	{
	  std::cout << "[HGCSimHitsAnalyzer][analyzeDigis] unable to find sector parameters for detId=0x" << hex << uint32_t(detId) << dec << " iSD=" << isd << " layer=" << layerKey << " sector=" << sector << std::endl;
	  continue;
	}
      
      //get local coordinates
      std::pair<float,float> xy = numberingSchemes_[isd]->getDDDConstants()->locateCell(cell,layer,subsector,true);
     
      //accumulate energy
      allSectors_[sectorKey][sector-1].digitize(rawDigi,xy.first,xy.second);
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

      int layer=detId.layer();
      int zpos=detId.zside();
      int layerKey(layer);
      if(zpos<0) layerKey *= -1;
      std::pair<int,int> sectorKey(isd,layerKey);
      int sector=detId.sector();
      int subsector=detId.subsector();
      int cell=detId.cell();
          
      if(allSectors_.find(sectorKey)==allSectors_.end()){
	std::cout << "[HGCSimHitsAnalyzer][analyzeDigis] unable to find layer parameters for detId=0x" << hex << uint32_t(detId) << dec << " iSD=" << isd << " layer=" << layerKey << std::endl;
	continue;
      }
      else if(allSectors_[sectorKey].size()<size_t(sector))
	{
	  std::cout << "[HGCSimHitsAnalyzer][analyzeDigis] unable to find sector parameters for detId=0x" << hex << uint32_t(detId) << dec << " iSD=" << isd << " layer=" << layerKey << " sector=" << sector << std::endl;
	  continue;
	}
      
      //get local coordinates
      std::pair<float,float> xy = numberingSchemes_[isd]->getDDDConstants()->locateCell(cell,layer,subsector,true);
     
      //accumulate energy
      allSectors_[sectorKey][sector-1].digitize(rawDigi,xy.first,xy.second);
    }
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSimHitsAnalyzer);
