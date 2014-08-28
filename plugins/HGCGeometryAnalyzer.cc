#include "UserCode/HGCanalysis/plugins/HGCGeometryAnalyzer.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Plane3D.h"
#include "CLHEP/Geometry/Vector3D.h"

#include "TF1.h"

#include <iostream>

#define SEP "\t\t\t\t"

using namespace std;

//
HGCGeometryAnalyzer::HGCGeometryAnalyzer( const edm::ParameterSet &iConfig )
  : testDone_(false)
{
  geometrySource_   = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");
  edm::Service<TFileService> fs;
  fs_=&fs;
}

//
HGCGeometryAnalyzer::~HGCGeometryAnalyzer()
{
}
  
//
void HGCGeometryAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  if(testDone_) return;

  //read geometry from event setup
  for(size_t i=0; i<geometrySource_.size(); i++)
    {
      edm::ESHandle<HGCalGeometry> hgcGeo;
      iSetup.get<IdealGeometryRecord>().get(geometrySource_[i],hgcGeo);

      const HGCalTopology &topo=hgcGeo->topology();
      const HGCalDDDConstants &dddConst=topo.dddConstants();

      int simLayers( dddConst.layers(false) );
      int recLayers( dddConst.layers(true) );

      //REPORT

      cout << geometrySource_[i] << " " << SEP   << "SIM "  <<  SEP  << "REC "  <<  SEP  << endl
	   << "-----------------------------------------------------------------------------" << endl
	   << "#layers " << SEP  << simLayers << SEP  << recLayers << SEP  << endl
	   << "-----------------------------------------------------------------------------" << endl
	   << "The following lines show the correspondence between sim and rec layers (#cells, #rows, cell size, b, t, h, alpha)" << endl
	   << "b/t=trapezoid half distance at the bottom/top  h=trapezoid height  alpha=half sector degree (0 if full sector)" << endl
	   << "Notice SIM/REC units are in mm/cm" << endl
	   << "-----------------------------------------------------------------------------" << endl;
      
      for(int ilay=1; ilay<=simLayers; ilay++)
	{
	  int simcells=dddConst.maxCells(ilay,false);
	  int simrows=dddConst.maxRows(ilay,false);
	  std::vector<HGCalDDDConstants::hgtrap>::const_iterator simModIt( dddConst.getFirstModule(false) );
	  for(int klay=1; klay<ilay; klay++) simModIt++;
	  std::cout.precision(4);
	  cout << " " << SEP  
	       << ilay << " : " 
	       << simcells << "/" 
	       << simrows << "/" 
	       << simModIt->cellSize << "/" 
	       << simModIt->bl << "/"
	       << simModIt->tl << "/"
	       << simModIt->h << "/"
	       << simModIt->alpha << SEP ;

	  std::pair<int,int>  simToReco=dddConst.simToReco(1,ilay,false);
	  std::vector<HGCalDDDConstants::hgtrap>::const_iterator recModIt( dddConst.getFirstModule(true) );
	  if(simToReco.second>0) 
	    {
	      int reccells=dddConst.maxCells(simToReco.second,true);
	      int recrows=dddConst.maxRows(simToReco.second,true);
	      for(int klay=1; klay<simToReco.second; klay++) recModIt++;
	      cout << simToReco.second << " : " 
		   << reccells << "/" 
		   << recrows << "/" 
		   << recModIt->cellSize << "/" 
		   << recModIt->bl << "/"
		   << recModIt->tl << "/"
		   << recModIt->h << "/"
		   << recModIt->alpha << SEP ;
	    }
	  else cout << "x" << SEP;
	  
	  cout << endl;

	  //test SIM->RECO coordinate assignment
	  if(simToReco.second<=0) continue;
	  if(ilay>1) continue;

	  TString name("layer"); name+= ilay; name +="_sd"; name+=i;

	  //sim
	  int ncellsy(2*simModIt->h/simModIt->cellSize);
	  int ncellsx(simModIt->tl/simModIt->cellSize);

	  TF1 *func=(*fs_)->make<TF1>("boundary_"+name,"[0]*x+[1]",0,ncellsx*simModIt->cellSize);
	  float func_slope     = (simModIt->alpha==0) ? (2*simModIt->h/(simModIt->tl-simModIt->bl)) : (simModIt->h/(simModIt->tl-simModIt->bl));
	  float func_offset    = -2*simModIt->h*simModIt->bl/(simModIt->tl-simModIt->bl)-simModIt->h;
	  func->SetParameter(0,func_slope);
	  func->SetParameter(1,func_offset);
	  func->SetLineWidth(2);
	  func->SetLineColor(1);

	  TH2F *simCellH = (*fs_)->make<TH2F>("simcell_"+ name, ";Local x [mm];Local y [mm];Sim cell number",
					      ncellsx,0,ncellsx*simModIt->cellSize,
					      ncellsy,-simModIt->h,ncellsy*simModIt->cellSize-simModIt->h);
	  TH2F *simixH = (TH2F *)simCellH->Clone("simix_"+name ); simixH = (*fs_)->make<TH2F>( *simixH );
	  TH2F *simiyH = (TH2F *)simCellH->Clone("simiy_"+name ); simiyH = (*fs_)->make<TH2F>( *simiyH );
	  
	  //reco 
	  TH2F *recCellH = (TH2F *)simCellH->Clone("reccell_"+name ); recCellH = (*fs_)->make<TH2F>( *recCellH );
	  TH2F *recixH = (TH2F *)simCellH->Clone("recix_"+name ); recixH = (*fs_)->make<TH2F>( *recixH );
	  TH2F *reciyH = (TH2F *)simCellH->Clone("reciy_"+name ); reciyH = (*fs_)->make<TH2F>( *reciyH );

	  //sim-reco diff
	  TH2F *dxH = (TH2F *)simCellH->Clone("dx_"+name ); dxH = (*fs_)->make<TH2F>( *dxH );
	  TH2F *dyH = (TH2F *)simCellH->Clone("dy_"+name ); dyH = (*fs_)->make<TH2F>( *dyH );
	  
	  for(int xbin=1; xbin<=simCellH->GetXaxis()->GetNbins(); xbin++)
	    for(int ybin=1; ybin<=simCellH->GetYaxis()->GetNbins(); ybin++)
	      {
		float x    = simCellH->GetXaxis()->GetBinCenter(xbin);
		float y    = simCellH->GetYaxis()->GetBinCenter(ybin);

		int simcell = dddConst.assignCell(x,y,ilay,1,false).second;
		simCellH->SetBinContent(xbin,ybin,simcell);
		int reccell = dddConst.assignCell(x/10.,y/10.,ilay,1,true).second;
		recCellH->SetBinContent(xbin,ybin,reccell);
		
		std::pair<int,int> ixy=dddConst.findCell(simcell,ilay,1,false);
		simixH->SetBinContent(xbin,ybin,ixy.first);
		simiyH->SetBinContent(xbin,ybin,ixy.second);
		
		std::pair<int,int> irecxy=dddConst.findCell(reccell,ilay,1,true);
		recixH->SetBinContent(xbin,ybin,irecxy.first);
		reciyH->SetBinContent(xbin,ybin,irecxy.second);
		
		//local coordinates (sim and reco)
		if(reccell>=0 && simcell>=0){
		  std::pair<float,float> simxy=dddConst.locateCell(simcell,ilay,1,false);
		  std::pair<float,float> recxy=dddConst.locateCell(reccell,ilay,1,true);
		  float dx=(simxy.first-recxy.first*10);
		  if(dx>10){
		    std::cout << dx << " " << simcell << ": "
			      << ixy.first << " " << ixy.second << " "
			      << "\t " << reccell << ": "
			      << irecxy.first << " " << irecxy.second << std::endl;
		  }
		  dxH->SetBinContent(xbin,ybin,dx);
		  float dy=(simxy.second-recxy.second*10);
		  dyH->SetBinContent(xbin,ybin,dy);
		}
	      }
	}
    }
  
  testDone_=true;
}




//define this as a plug-in
DEFINE_FWK_MODULE(HGCGeometryAnalyzer);
