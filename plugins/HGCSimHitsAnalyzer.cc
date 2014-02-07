#include "UserCode/HGCanalysis/plugins/HGCSimHitsAnalyzer.h"

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


#include <iostream>

using namespace std;

//
HGCSimHitsAnalyzer::HGCSimHitsAnalyzer( const edm::ParameterSet &iConfig ) : geometryDefined_(false), numberingScheme_(0)
{
  ddViewName_  = iConfig.getUntrackedParameter<std::string>("ddViewName", "");
  eeHits_      = iConfig.getUntrackedParameter<std::string>("eeHits",     "HGCHitsEE");
  heHits_      = iConfig.getUntrackedParameter<std::string>("heHits",     "HGCHitsHE");
  genSource_   = iConfig.getUntrackedParameter<std::string>("genSource",  "genParticles");

  edm::Service<TFileService> fs;
  t_=fs->make<TTree>("HGC","Event Summary");
  t_->Branch("run",      &simEvt_.run,      "run/I");
  t_->Branch("lumi",     &simEvt_.lumi,     "lumi/I");
  t_->Branch("event",    &simEvt_.event,    "event/I");
  t_->Branch("nee",      &simEvt_.nee,      "nee/I");
  t_->Branch("ee_zp",     simEvt_.ee_zp, "ee_zp[nee]/I");
  t_->Branch("ee_layer",  simEvt_.ee_layer, "ee_layer[nee]/I");
  t_->Branch("ee_module", simEvt_.ee_module, "ee_module[nee]/I");
  t_->Branch("ee_subsec", simEvt_.ee_subsec, "ee_subsec[nee]/I");
  t_->Branch("ee_cell",   simEvt_.ee_cell,   "ee_cell[nee]/I");
  t_->Branch("ee_edep",   simEvt_.ee_edep,  "ee_edep[nee]/F");
  t_->Branch("ee_t",      simEvt_.ee_t,     "ee_t[nee]/F");
  t_->Branch("ee_x",      simEvt_.ee_x,     "ee_x[nee]/F");
  t_->Branch("ee_y",      simEvt_.ee_y,     "ee_y[nee]/F");
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

  //get hits and gen level information
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  iEvent.getByLabel(edm::InputTag(genSource_), genParticles);

  edm::Handle<edm::PCaloHitContainer> caloHitsEE;
  iEvent.getByLabel(edm::InputTag("g4SimHits",eeHits_),caloHitsEE); 

  edm::Handle<edm::PCaloHitContainer> caloHitsHE;
  iEvent.getByLabel(edm::InputTag("g4SimHits",heHits_),caloHitsHE); 

  simEvt_.run    = iEvent.id().run();
  simEvt_.lumi   = iEvent.luminosityBlock();
  simEvt_.event  = iEvent.id().event();

  simEvt_.nee=0;
  analyzeEEHits(caloHitsEE,genParticles);
  
  if(simEvt_.nee>0)  t_->Fill();
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
  numberingScheme_=new HGCNumberingScheme(gpar);

  
  //parse the DD for sensitive volumes
  DDExpandedView eview(cview);
  std::map<DDExpandedView::nav_type,int> idMap;
  do {
    const DDLogicalPart &logPart=eview.logicalPart();
    std::string name=logPart.name();

    //only EE sensitive volumes for the moment
    if(name.find("Sensitive")==std::string::npos) continue;
    if(name.find("HGCalEE")==std::string::npos) continue;
    
    size_t pos=name.find("Sensitive")+9;
    int layer=atoi(name.substr(pos,name.size()).c_str());
   
    //save half height and widths for the trapezoid
    if(eeSVpars_.find(layer)!=eeSVpars_.end()) continue; 
    std::vector<double> solidPars=eview.logicalPart().solid().parameters();
    std::cout << logPart << endl;
    for(size_t i=0; i<solidPars.size(); i++) std::cout << solidPars[i] << " ";
    std::cout << endl;

    std::vector<double> layerPars;
    layerPars.push_back( solidPars[3] ); //height
    layerPars.push_back( solidPars[4] ); //bottom
    layerPars.push_back( solidPars[5] ); //top
    eeSVpars_[ layer ] = layerPars;
  }while(eview.next() );

  //all done here
  return true;
}

//
void HGCSimHitsAnalyzer::analyzeEEHits(edm::Handle<edm::PCaloHitContainer> &caloHits, edm::Handle<edm::View<reco::Candidate> > &gen)
{
  if(!caloHits.isValid()) return;
  for(edm::PCaloHitContainer::const_iterator hit_it = caloHits->begin(); hit_it != caloHits->end(); ++hit_it) 
    {
      HGCEEDetId detId(hit_it->id());

      int layer=detId.layer();
      if(eeSVpars_.find(layer) == eeSVpars_.end()){
	std::cout << "[HGCSimHitsAnalyzer][analyzeEEHits] unable to find layer parameters for detId=0x" << hex << uint32_t(detId) << dec << std::endl;
	//continue;
      }
      
      int cell=detId.cell();
      std::pair<float,float> xy = numberingScheme_->getLocalCoords(cell,
								   numberingScheme_->getCellSize(),
								   eeSVpars_[layer][HALF_H],
								   eeSVpars_[layer][HALF_B],
								   eeSVpars_[layer][HALF_T]);
      int subsector=detId.subsector();
      if(subsector==0) xy.first *=-1;

      int module=detId.module();

      simEvt_.ee_zp[simEvt_.nee]     = detId.zside();
      simEvt_.ee_layer[simEvt_.nee]  = layer;
      simEvt_.ee_module[simEvt_.nee] = module;
      simEvt_.ee_subsec[simEvt_.nee]   = detId.subsector();
      simEvt_.ee_cell[simEvt_.nee]   = detId.cell();
      simEvt_.ee_edep[simEvt_.nee]   = hit_it->energy();
      simEvt_.ee_t[simEvt_.nee]      = hit_it->time();
      simEvt_.ee_x[simEvt_.nee]      = xy.first;
      simEvt_.ee_y[simEvt_.nee]      = xy.second;
      simEvt_.nee++;

//       std::cout << hex << "0x" << uint32_t(detId) << dec 
// 		<< " | isFwd=" << detId.isForward()
// 		<< " isEE=" << detId.isEE() 
// 		<< " zside=" << detId.zside() 
// 		<< " layer=" << detId.layer() 
// 		<< " subsec=" << detId.subsector() 
// 		<< " module=" << detId.module()
// 		<< " cell=" << detId.cell() 
// 		<< " | den=" << den
// 		<< std::endl;
    }
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSimHitsAnalyzer);
