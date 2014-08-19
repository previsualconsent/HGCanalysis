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
HGCSimHitsAnalyzer::HGCSimHitsAnalyzer( const edm::ParameterSet &iConfig ) : 
  saveGenParticles_(true),
  saveG4_(true),
  saveTkExtrapol_(true),
  piTkPropagator_(0)
{
  //configure analyzer
  saveGenParticles_ = iConfig.getUntrackedParameter< bool > ("saveGenParticles"); 
  genSource_        = iConfig.getUntrackedParameter<std::string>("genSource");
  saveG4_           = iConfig.getUntrackedParameter< bool > ("saveG4");
  g4TracksSource_   = iConfig.getUntrackedParameter<std::string>("g4TracksSource");
  g4VerticesSource_ = iConfig.getUntrackedParameter<std::string>("g4VerticesSource");
  saveTkExtrapol_   = iConfig.getUntrackedParameter< bool > ("saveTkExtrapol");
  trackSource_      = iConfig.getUntrackedParameter<std::string>("trackSource");
  geometrySource_   = iConfig.getParameter< std::vector<std::string> >("geometrySource");
  hitCollections_   = iConfig.getUntrackedParameter< std::vector<std::string> >("hitCollections");
  digiCollections_  = iConfig.getUntrackedParameter< std::vector<std::string> >("digiCollections");

  //init tree
  edm::Service<TFileService> fs;
  t_=fs->make<TTree>("HGC","Event Summary");
  initHGCSimulationEventTree(t_,simEvt_);
}

//
HGCSimHitsAnalyzer::~HGCSimHitsAnalyzer()
{
}

//
void HGCSimHitsAnalyzer::analyzeGenParticles(edm::Handle<edm::View<reco::Candidate> > &genParticles)
{
  simEvt_.ngen=0;
  if(!genParticles.isValid()) return;
  if(genParticles->size()==0) return;

  for(size_t i = 0; i < genParticles->size(); ++ i)
    {
      const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[i] );
      simEvt_.gen_id[simEvt_.ngen]=p.pdgId();
      simEvt_.gen_pt[simEvt_.ngen]=p.pt();
      simEvt_.gen_eta[simEvt_.ngen]=p.eta();
      simEvt_.gen_phi[simEvt_.ngen]=p.phi();
      simEvt_.gen_en[simEvt_.ngen]=p.energy();
      simEvt_.ngen++;
      if(simEvt_.ngen>=MAXGENPEREVENT) break;
    }
}

//
void HGCSimHitsAnalyzer::analyzeG4information(edm::Handle<edm::SimTrackContainer> &SimTk,edm::Handle<edm::SimVertexContainer> &SimVtx)
{
  simEvt_.ng4=0;
  if(!SimTk.isValid() || !SimVtx.isValid()) return;
  if(SimTk->size()==0) return;
  
  //the primary
  const SimTrack &primtk=SimTk->at(0);
  const math::XYZTLorentzVectorD& primpos = SimVtx->at(0).position();
  const math::XYZTLorentzVectorD &primp4  = primtk.momentum() ;
  simEvt_.g4_id [0] = primtk.type();
  simEvt_.g4_vtx[0] = primpos.x();
  simEvt_.g4_vty[0] = primpos.y();
  simEvt_.g4_vtz[0] = primpos.z();
  simEvt_.g4_en [0] = primp4.energy();
  simEvt_.g4_eta[0] = primp4.eta();
  simEvt_.g4_phi[0] = primp4.phi();
  simEvt_.ng4++;
  
  //compute energy loss in tracker 
  // save secondaries in the calorimeter
  simEvt_.g4_dEnInTracker=0;
  for (unsigned int isimvtx = 1; isimvtx < SimVtx->size() ; isimvtx++ )
    {
      const math::XYZTLorentzVectorD& pos=SimVtx->at(isimvtx).position();
      bool isInTracker( fabs(pos.z())< 320); //yes...this is hardcoded

      for (unsigned int isimtk = 0; isimtk < SimTk->size() ; isimtk++ )
	{
	  //get all tracks for this vertex
	  const SimTrack &tk=SimTk->at(isimtk);
	  int vtxIdx=tk.vertIndex();
	  if(vtxIdx<0) continue;
	  if(vtxIdx!=(int) isimvtx) continue;

	  //if primary and in the tracker update its energy
	  if( isInTracker && tk.genpartIndex()==primtk.genpartIndex() ) 
	    {
	      simEvt_.g4_dEnInTracker = tk.momentum().energy();
	    }

	  //if in the calorimeter save everything
	  else if(!isInTracker && simEvt_.ng4<MAXG4PEREVENT)
	    {
	      const math::XYZTLorentzVectorD &p4 = tk.momentum() ;
	      simEvt_.g4_id[simEvt_.ngen] = tk.type();
	      simEvt_.g4_vtx[simEvt_.ngen] = pos.x();
	      simEvt_.g4_vty[simEvt_.ngen] = pos.y();
	      simEvt_.g4_vtz[simEvt_.ngen] = pos.z();
	      simEvt_.g4_en[simEvt_.ngen] = p4.energy();
	      simEvt_.g4_eta[simEvt_.ngen] = p4.eta();
	      simEvt_.g4_phi[simEvt_.ngen] = p4.phi();
	      simEvt_.ngen++;
	    }
	}
    }

  //subtract incoming energy
  simEvt_.g4_dEnInTracker=-simEvt_.g4_en[0];
}

//
void HGCSimHitsAnalyzer::analyzeTrackingInformation(edm::Handle<reco::TrackCollection> &tracks,
						    edm::ESHandle<TrackerGeometry> &tkGeom,
						    edm::ESHandle<MagneticField> &bField,
						    std::map<int,const HGCalGeometry *> &hgcGeometries)
{
  simEvt_.ntk=0;
  if(!tracks.isValid() || !tkGeom.isValid() || !bField.isValid()) return;
  if(tracks->size()==0) return;
  
  //init track propagator and surfaces to be used, if not yet done
  if(piTkPropagator_==0)
    {
      cout << "[HGCSimHitsAnalyzer][analyzeTrackingInformation] starting tracking propagators" << endl;
      piTkPropagator_ = new PropagatorWithMaterial(alongMomentum,0.1396, bField.product());

      Surface::RotationType rot; //unit rotation matrix
      for(std::map<int,const HGCalGeometry *>::iterator it = hgcGeometries.begin(); it!= hgcGeometries.end(); it++)
	{
	  std::cout << "\t HGC subdet: " << it->first << std::endl;
	  const HGCalDDDConstants &dddCons=it->second->topology().dddConstants();
	  std::map<float,float> zrhoCoord;
	  std::vector< ReferenceCountingPointer<BoundDisk> > iMinusSurfaces, iPlusSurfaces;
	  std::vector<HGCalDDDConstants::hgtrform>::const_iterator firstLayerIt = dddCons.getFirstTrForm();
	  std::vector<HGCalDDDConstants::hgtrform>::const_iterator lastLayerIt = dddCons.getLastTrForm();
	  for(std::vector<HGCalDDDConstants::hgtrform>::const_iterator layerIt=firstLayerIt; layerIt!=lastLayerIt; layerIt++)
	    {
	      float Z(fabs(layerIt->h3v.z()));
	      float Radius(dddCons.getLastModule(true)->tl+layerIt->h3v.perp());
	      zrhoCoord[Z]=Radius;
	    }
	  for(std::map<float,float>::iterator it=zrhoCoord.begin(); it!= zrhoCoord.end(); it++)
	    {
	      float Z(it->first);
	      float Radius(it->second);
	      std::cout << "\t z=" << Z << std::flush;
	      iMinusSurfaces.push_back(ReferenceCountingPointer<BoundDisk> ( new BoundDisk( Surface::PositionType(0,0,-Z), rot, new SimpleDiskBounds( 0, Radius, -0.001, 0.001))));
	      iPlusSurfaces.push_back(ReferenceCountingPointer<BoundDisk> ( new BoundDisk( Surface::PositionType(0,0,+Z), rot, new SimpleDiskBounds( 0, Radius, -0.001, 0.001))));
	    }
	  minusSurface_[it->first] = iMinusSurfaces;
	  plusSurface_[it->first] = iPlusSurfaces;
	  std::cout << " | total " << minusSurface_[it->first].size() << " layers for subdet #" << it->first << std::endl;
	}
    }


  //analyze tracks
  for(std::vector<reco::Track>::const_iterator tIt = tracks->begin(); tIt != tracks->end(); tIt++)
    {
      //it should be safe to cut away tracks which are not pointing to the endcaps
      if(fabs(tIt->eta())<1.45 || fabs(tIt->eta())>3.5) continue;
      simEvt_.tk_pt[simEvt_.ntk] = tIt->pt();
      simEvt_.tk_eta[simEvt_.ntk] = tIt->eta();
      simEvt_.tk_phi[simEvt_.ntk] = tIt->phi();
      simEvt_.tk_chi2[simEvt_.ntk] = tIt->normalizedChi2();
      simEvt_.tk_nhits[simEvt_.ntk] = tIt->hitPattern().trackerLayersWithMeasurement();

      const TrajectoryStateOnSurface myTSOS = trajectoryStateTransform::outerStateOnSurface(*tIt, *(tkGeom.product()), bField.product());
      std::map<int, std::vector<ReferenceCountingPointer<BoundDisk> > >::iterator itBegin( myTSOS.globalPosition().z()>0 ? plusSurface_.begin() : minusSurface_.begin() );
      std::map<int, std::vector<ReferenceCountingPointer<BoundDisk> > >::iterator itEnd  ( myTSOS.globalPosition().z()>0 ? plusSurface_.end()   : minusSurface_.end() );
      int iextrapol(0);
      for(std::map<int, std::vector<ReferenceCountingPointer<BoundDisk> > >::iterator it = itBegin; it!=itEnd; it++)
	{
	  for(size_t ilayer=0; ilayer<it->second.size(); ilayer++)
	    {
	      TrajectoryStateOnSurface piStateAtSurface = piTkPropagator_->propagate (myTSOS, *((it->second)[ilayer]) );
	      if(piStateAtSurface.isValid())
		{
		  GlobalPoint pt=piStateAtSurface.globalPosition();
		  simEvt_.tk_extrapol_sd[simEvt_.ntk][iextrapol]=it->first;
		  simEvt_.tk_extrapol_layer[simEvt_.ntk][iextrapol]=ilayer;
		  simEvt_.tk_extrapol_x[simEvt_.ntk][iextrapol]=pt.x();
		  simEvt_.tk_extrapol_y[simEvt_.ntk][iextrapol]=pt.y();
		  iextrapol++;
		}
	    }
	}
      simEvt_.ntk++;
      if(simEvt_.ntk>=MAXTKSPEREVENT) break;
    }
}

  
//
void HGCSimHitsAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  cout << "[HGCSimHitsAnalyzer::analyze] start" << endl;
  
  //event header
  simEvt_.run    = iEvent.id().run();
  simEvt_.lumi   = iEvent.luminosityBlock();
  simEvt_.event  = t_->GetEntriesFast()+1;
  
  //generator level particles
  if(saveGenParticles_)
    {
      edm::Handle<edm::View<reco::Candidate> > genParticles;
      iEvent.getByLabel(edm::InputTag(genSource_), genParticles);
      analyzeGenParticles(genParticles);
    }
  
  //Geant4 information
  if(saveG4_)
    {
      edm::Handle<edm::SimTrackContainer> SimTk;
      iEvent.getByLabel(g4TracksSource_,SimTk);
      edm::Handle<edm::SimVertexContainer> SimVtx;
      iEvent.getByLabel(g4VerticesSource_,SimVtx);
      analyzeG4information(SimTk,SimVtx);
    }
  
  //read geometry from event setup
  std::map<int,const HGCalGeometry *> hgcGeometries;
  for(size_t i=0; i<geometrySource_.size(); i++)
    {
      edm::ESHandle<HGCalGeometry> hgcGeo;
      iSetup.get<IdealGeometryRecord>().get(geometrySource_[i],hgcGeo);
      hgcGeometries[i]=hgcGeo.product();
    }

  //save tracker extrapolation
  if(saveTkExtrapol_)
    {
      edm::ESHandle<MagneticField> bField;
      iSetup.get<IdealMagneticFieldRecord>().get(bField);
      edm::Handle<reco::TrackCollection> tracks;
      iEvent.getByLabel(edm::InputTag(trackSource_), tracks);
      edm::ESHandle<TrackerGeometry> tkGeom;
      iSetup.get<TrackerDigiGeometryRecord>().get( tkGeom );
      analyzeTrackingInformation(tracks,tkGeom,bField,hgcGeometries);
    }


  //hits, digis
  simEvt_.nhits=0;
  for(size_t i=0; i<hitCollections_.size(); i++)
    {
      edm::Handle<edm::PCaloHitContainer> caloHits;
      iEvent.getByLabel(edm::InputTag("g4SimHits",hitCollections_[i]),caloHits); 
      analyzeHits(i,caloHits,hgcGeometries[i]);
      
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
      for(HGCSimHitDataAccumulator::iterator dit=simHitData_[isd].begin(); dit!=simHitData_[isd].end(); dit++)
	{
	  Float_t edep(0);
	  for(size_t itime=0; itime<=8; itime++) edep+=dit->second[itime]*1e6;
	  if(edep<=0) continue;
	  if(simEvt_.nhits>=MAXHGCHITSPEREVENT) break;

	  uint32_t rawDetId = dit->first;
	  uint32_t layer    = (rawDetId>>19) & 0x1F;	  
	  simEvt_.hit_type [simEvt_.nhits]=isd;
	  simEvt_.hit_layer[simEvt_.nhits]=layer;
	  simEvt_.hit_edep [simEvt_.nhits]=edep;
	  for(size_t itime=0; itime<=8; itime++) simEvt_.hit_edep_sample[simEvt_.nhits][itime]=dit->second[itime]*1e6;
	  UShort_t adc=(UShort_t)dit->second[9];
	  simEvt_.digi_adc [simEvt_.nhits]=adc;
	  
	  //get position
	  const GlobalPoint pos( std::move( hgcGeometries[isd]->getPosition( rawDetId) ) );
 	  simEvt_.hit_x[simEvt_.nhits]=pos.x();
 	  simEvt_.hit_y[simEvt_.nhits]=pos.y();
 	  simEvt_.hit_z[simEvt_.nhits]=pos.z();
 	  simEvt_.hit_eta[simEvt_.nhits]=pos.eta();
 	  simEvt_.hit_phi[simEvt_.nhits]=pos.phi();
	
	  simEvt_.nhits++;
	}
    }
 
  //fill tree
  if(simEvt_.nhits>0)  t_->Fill();
}



//
void HGCSimHitsAnalyzer::analyzeHits(size_t isd,edm::Handle<edm::PCaloHitContainer> &caloHits,const HGCalGeometry *geom)
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

      //gang SIM->RECO cells
      int layer(simId.layer()), cell(simId.cell());
      float zPos(0.);
      const HGCalTopology &topo=geom->topology();
      const HGCalDDDConstants &dddConst=topo.dddConstants();
      zPos=dddConst.getFirstTrForm()->h3v.z();
      std::pair<int,int> recoLayerCell=dddConst.simToReco(cell,layer,topo.detectorType());
      cell  = recoLayerCell.first;
      layer = recoLayerCell.second;
      if(layer<0) continue;
       
      //assign the RECO DetId
      uint32_t id( (isd==0) ?
		   (uint32_t)HGCEEDetId(ForwardSubdetector(mySubDet),simId.zside(),layer,simId.sector(),simId.subsector(),cell):
		   (uint32_t)HGCHEDetId(ForwardSubdetector(mySubDet),simId.zside(),layer,simId.sector(),simId.subsector(),cell) );
      
      //hit time: [time()]=ns  [zPos]=cm [CLHEP::c_light]=mm/ns
      //for now accumulate in buckets of 6.25ns = 4 time samples each 25 ns
      int itime=floor( (hit_it->time()-zPos/(0.1*CLHEP::c_light))/6.25);
      
      HGCSimHitDataAccumulator::iterator simHitIt=simHitData_[isd].find(id);
      if(simHitIt==simHitData_[isd].end())
        {
          HGCSimHitData baseData;
	  baseData.fill(0.);
          simHitData_[isd][id]=baseData;
          simHitIt=simHitData_[isd].find(id);
	}
      if( itime > 8) continue; //last will be reserved for the ADC counts
      (simHitIt->second)[itime] += hit_it->energy();
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
          HGCSimHitData baseData;
	  baseData.fill(0.);
          simHitData_[isd][detId]=baseData;
          simHitIt=simHitData_[isd].find(detId);
        }
      (simHitIt->second)[9] = rawDigi;
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
          HGCSimHitData baseData;
	  baseData.fill(0.);
          simHitData_[isd][detId]=baseData;
          simHitIt=simHitData_[isd].find(detId);
        }
      (simHitIt->second)[9] = rawDigi;
    }
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSimHitsAnalyzer);
