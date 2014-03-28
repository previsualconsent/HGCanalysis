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

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"

#include "TVector2.h"

#include <iostream>

using namespace std;

//
HGCSimHitsAnalyzer::HGCSimHitsAnalyzer( const edm::ParameterSet &iConfig ) : geometryDefined_(false)
{
  ddViewName_     = iConfig.getUntrackedParameter<std::string>("ddViewName", "");
  genSource_      = iConfig.getUntrackedParameter<std::string>("genSource",  "genParticles");
  hitCollections_ = iConfig.getUntrackedParameter< std::vector<std::string> >("hitCollections");
  cellSizePars_   = iConfig.getUntrackedParameter< std::vector<int> >("cellSizePars");
  sdTags_         = iConfig.getUntrackedParameter< std::vector<std::string> >("sdTags");

  edm::Service<TFileService> fs;
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
  t_->Branch("hit_zp",     simEvt_.hit_zp,     "hit_zp[nhits]/I");
  t_->Branch("hit_layer",  simEvt_.hit_layer,  "hit_layer[nhits]/I");
  t_->Branch("hit_sec",    simEvt_.hit_sec,    "hit_sec[nhits]/I");
  t_->Branch("hit_subsec", simEvt_.hit_subsec, "hit_subsec[nhits]/I");
  t_->Branch("hit_cell",   simEvt_.hit_cell,   "hit_cell[nhits]/I");
  t_->Branch("hit_edep",   simEvt_.hit_edep,   "hit_edep[nhits]/F");
  t_->Branch("hit_t",      simEvt_.hit_t,      "hit_t[nhits]/F");
  t_->Branch("hit_x",      simEvt_.hit_x,      "hit_x[nhits]/F");
  t_->Branch("hit_y",      simEvt_.hit_y,      "hit_y[nhits]/F");
  t_->Branch("hit_gx",     simEvt_.hit_gx,     "hit_gx[nhits]/F");
  t_->Branch("hit_gy",     simEvt_.hit_gy,     "hit_gy[nhits]/F");
  t_->Branch("hit_gz",     simEvt_.hit_gz,     "hit_gz[nhits]/F");
  
  for(size_t i=0; i<hitCollections_.size(); i++)
    { 
      TString prefix("sens"); prefix+=i;
      sensNameH_   .push_back( fs->make<TH1F>(prefix+"_name",            ";Sensitive detector;Present",1,0,1) );
      sensNameH_[i]->GetXaxis()->SetBinLabel( 1, hitCollections_[i].c_str() );
      sensHeightH_ .push_back( fs->make<TH1F>(prefix+"_HalfHeight",      ";Layer;Half height [mm]",       400,-200.,200.) );
      sensBottomH_ .push_back( fs->make<TH1F>(prefix+"_BottomHalfWidth", ";Layer;Bottom half width [mm]", 400,-200.,200.) );
      sensTopH_    .push_back( fs->make<TH1F>(prefix+"_TopHalfWidth",    ";Layer;Top half width [mm]",    400,-200.,200.) );
      sensBasePhiH_.push_back( fs->make<TH1F>(prefix+"_BasePhi",         ";Layer;Baseline #phi [rad]",    400,-200.,200.) );
      sensTranslXH_.push_back( fs->make<TH1F>(prefix+"_TranslX",         ";Layer;Translation X [mm]",     400,-200.,200.) );
      sensTranslYH_.push_back( fs->make<TH1F>(prefix+"_TranslY",         ";Layer;Translation Y [mm]",     400,-200.,200.) );
      sensTranslZH_.push_back( fs->make<TH1F>(prefix+"_TranslZ",         ";Layer;Translation Z [mm]",     400,-200.,200.) );
    }

  sensSVpars_.resize( hitCollections_.size() );
}

//
HGCSimHitsAnalyzer::~HGCSimHitsAnalyzer()
{
}

//
void HGCSimHitsAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  cout << "[HGCSimHitsAnalyzer::analyze] start" << endl;

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

  //get gen level information
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  iEvent.getByLabel(edm::InputTag(genSource_), genParticles);
  simEvt_.ngen=0;
  for(size_t i = 0; i < genParticles->size(); ++ i)
    {
      const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticles)[i] );
      if(p.status()!=1) continue;
      if(abs(p.pdgId())!=11 && abs(p.pdgId())!=13) continue;
      simEvt_.gen_id[simEvt_.ngen]=p.pdgId();
      simEvt_.gen_pt[simEvt_.ngen]=p.pt();
      simEvt_.gen_eta[simEvt_.ngen]=p.eta();
      simEvt_.gen_phi[simEvt_.ngen]=p.phi();
      simEvt_.gen_en[simEvt_.ngen]=p.energy();
      simEvt_.ngen++;
    }

  simEvt_.run    = iEvent.id().run();
  simEvt_.lumi   = iEvent.luminosityBlock();
  simEvt_.event  = t_->GetEntriesFast()+1; //iEvent.id().event();

  simEvt_.nhits=0;
  for(size_t i=0; i<hitCollections_.size(); i++)
    {
      edm::Handle<edm::PCaloHitContainer> caloHits;
      iEvent.getByLabel(edm::InputTag("g4SimHits",hitCollections_[i]),caloHits); 
      analyzeHits(i,caloHits,genParticles);
    }
  
  if(simEvt_.nhits>0)  t_->Fill();
}

//
bool HGCSimHitsAnalyzer::defineGeometry(edm::ESTransientHandle<DDCompactView> &ddViewH)
{
  if(!ddViewH.isValid()) {
    std::cout << "[HGCSimHitsAnalyzer][defineGeometry] invalid DD handle" << std::endl;
    return false;
  }
  
  const DDCompactView &cview=*ddViewH;

  //get geometry parameters from DDD (cell size, layer limits, etc.)
  DDSpecificsFilter filter0;
  DDValue ddv0("Volume", "HGC", 0);
  filter0.setCriteria(ddv0, DDSpecificsFilter::equals);
  DDFilteredView fv0(cview);
  fv0.addFilter(filter0);
  fv0.firstChild();
  DDsvalues_type sv0(fv0.mergedSpecifics());
  std::vector<double> gpar;
  DDValue ddv1("GeomParHGC");
  if(DDfetch(&sv0,ddv1))  gpar = ddv1.doubles();
  if(gpar.size()==0){
    std::cout << "[HGCSimHitsAnalyzer][defineGeometry] failed to retrieve geometry parameters from DDD" << std::endl;
    return false;
  }

  //instantiate the numbering schemes
  for(size_t i=0; i<cellSizePars_.size(); i++)
    {
      std::vector<double> igpar(1, gpar[ cellSizePars_[i] ] );
      numberingSchemes_.push_back( new HGCNumberingScheme(igpar) );
      std::cout << "[HGCSimHitsAnalyzer][defineGeometry] cell size for " << i << "-th hit collection is " << igpar[0] << std::endl;
    }
  
  //parse the DD for sensitive volumes
  DDExpandedView eview(cview);
  std::map<DDExpandedView::nav_type,int> idMap;
  do {
    const DDLogicalPart &logPart=eview.logicalPart();
    std::string name=logPart.name();

    //only EE sensitive volumes for the moment
    size_t isd(999999);
    for(size_t i=0; i<sdTags_.size(); i++)
      {
	if(name.find(sdTags_[i])==std::string::npos) continue;
	isd=i;
	break;
      }
    if(isd==999999) continue; // a bit ugly ...


    //this is common
    size_t pos=name.find("Sensitive")+9;
    int layer=atoi(name.substr(pos,name.size()).c_str());

    //translation and rotation for this part
    //note that the inverse rotation matrix elements are:
    //[xrot^T yrot^T zrot^T]
    //cf. http://root.cern.ch/root/html/src/ROOT__Math__Rotation3D.h.html#twSZND
    DDTranslation transl=eview.translation();
    DDRotationMatrix rot=eview.rotation();
 
    //set key to -1 if in negative z axis
    int layerKey(layer);
    if( transl.z()<0 ) layerKey *= -1;
   
//     DDExpandedNode parent = eview.geoHistory()[ eview.geoHistory().size()-2 ];
//     DDTranslation parentTransl = parent.absTranslation();
//     DDRotationMatrix parentRot = parent.absRotation();
//     transl = parentRot.Inverse()*(transl - parentTransl );
//     rot = parentRot.Inverse()*rot;
//     rot = rot.Inverse(); 

    DD3Vector xrot, yrot, zrot;
    rot.GetComponents(xrot,yrot,zrot);
    double basePhi=TMath::ATan2(-yrot.x(),xrot.x());
 
    std::vector<double> solidPars=eview.logicalPart().solid().parameters();
    
    SectorGeometry_t isector;
    isector.basePhi=basePhi;
    isector.halfHeight=solidPars[3];
    isector.halfBottom=solidPars[4];
    isector.halfTop=solidPars[5];
    isector.halfWidth=solidPars[0];
    isector.globalX=transl.x();
    isector.globalY=transl.y();
    isector.globalZ=transl.z();
    isector.xx=xrot.x(); isector.xy=yrot.x(); isector.xz=zrot.x();
    isector.yx=xrot.y(); isector.yy=yrot.y(); isector.yz=zrot.y();
    isector.zx=xrot.z(); isector.zy=yrot.z(); isector.zz=zrot.z();

    //save half height and widths for the trapezoid
    if(sensSVpars_[isd].find(layerKey)==sensSVpars_[isd].end()) 
      {
	std::vector<SectorGeometry_t> sectorPars(1,isector);
	sensSVpars_[isd][layerKey]=sectorPars;
	sensNameH_[isd]      ->Fill(0);
 	sensHeightH_[isd]    ->Fill(layerKey, isector.halfHeight);
 	sensBottomH_[isd]    ->Fill(layerKey, isector.halfBottom );
 	sensTopH_[isd]       ->Fill(layerKey, isector.halfTop );
 	sensTranslXH_[isd]   ->Fill(layerKey, isector.globalX);
 	sensTranslYH_[isd]   ->Fill(layerKey, isector.globalY);
 	sensTranslZH_[isd]   ->Fill(layerKey, isector.globalZ);
 	sensBasePhiH_[isd]   ->Fill(layerKey, isector.basePhi);
      }
    else
      sensSVpars_[isd][layerKey].push_back(isector);
    
  }while(eview.next() );

  //all done here
  return true;
}

//
void HGCSimHitsAnalyzer::analyzeHits(size_t isd,edm::Handle<edm::PCaloHitContainer> &caloHits, edm::Handle<edm::View<reco::Candidate> > &gen)
{
  //check inputs
  if(!caloHits.isValid()) return;
  if(sensSVpars_.size()<=isd)
    {
      std::cout << "[HGCSimHitsAnalyzer][analyzeHits] unable to find SD parameters for iSD=" << isd << std::endl;
      return;
    }
  
  //analyze hits
  for(edm::PCaloHitContainer::const_iterator hit_it = caloHits->begin(); hit_it != caloHits->end(); ++hit_it) 
    {
      HGCEEDetId detId(hit_it->id());
      
      int layer=detId.layer();
      int zpos=detId.zside();

      //check layer
      int layerKey(layer);
      if(zpos<0) layerKey *= -1;
      if(sensSVpars_[isd].find(layerKey) == sensSVpars_[isd].end()){
	std::cout << "[HGCSimHitsAnalyzer][analyzeHits] unable to find layer parameters for detId=0x" << hex << uint32_t(detId) << dec << " iSD=" << isd << " layer=" << layerKey << std::endl;
	continue;
      }
      
      //check sector
      int sector=detId.sector();
      if(sensSVpars_[isd][layerKey].size()<size_t(sector))
	{
	  std::cout << "[HGCSimHitsAnalyzer][analyzeHits] unable to find sector parameters for detId=0x" << hex << uint32_t(detId) << dec << " iSD=" << isd << " layer=" << layerKey << " sector=" << sector << std::endl;
	  continue;
	}
      
      //get the geometry info
      SectorGeometry_t &isector=sensSVpars_[isd][layerKey][sector-1];

      int cell=detId.cell();
      std::pair<float,float> xy = numberingSchemes_[isd]->getLocalCoords(cell,
									 numberingSchemes_[isd]->getCellSize(),
									 isector.halfHeight,
									 isector.halfBottom,
									 isector.halfTop);
      float localX(xy.first), localY(xy.second);
      int subsector=detId.subsector();
      localX *= subsector;
      
      //hit type
      simEvt_.hit_type[simEvt_.nhits]     = isd;

      //det id info
      simEvt_.hit_zp[simEvt_.nhits]     = zpos;
      simEvt_.hit_layer[simEvt_.nhits]  = layer;
      simEvt_.hit_sec[simEvt_.nhits]    = sector;
      simEvt_.hit_subsec[simEvt_.nhits] = detId.subsector();
      simEvt_.hit_cell[simEvt_.nhits]   = detId.cell();

      //energy
      simEvt_.hit_edep[simEvt_.nhits]   = hit_it->energy();

      //time
      simEvt_.hit_t[simEvt_.nhits]      = hit_it->time();

      //local coordinates (relative to the center of the sensitive detector)
      simEvt_.hit_x[simEvt_.nhits]      = localX;
      simEvt_.hit_y[simEvt_.nhits]      = localY;

      //global coordinates 
      float rho=sqrt(pow(isector.globalX,2)+pow(isector.globalY,2));
      float rotLocalX=localX;
      float rotLocalY=localY+rho;
      simEvt_.hit_gx[simEvt_.nhits]     = rotLocalX*isector.xx+rotLocalY*isector.xy;
      simEvt_.hit_gy[simEvt_.nhits]     = rotLocalX*isector.yx+rotLocalY*isector.yy;
      simEvt_.hit_gz[simEvt_.nhits]     = isector.globalZ;

      //increment array
      simEvt_.nhits++;
    }
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSimHitsAnalyzer);
