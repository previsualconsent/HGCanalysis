{
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile *_file0 = TFile::Open("/tmp/psilva/HGCGeometry.root");
  TString outDir("~/www/HGCal/Geometry/v4/");

  TCanvas *c=new TCanvas("c","c",500,500);
  c->SetRightMargin(0.2);
  TString dists[]={"simcell","reccell","dx","dy","recix","reciy","simix","simiy"};
  for(int ilayer=1; ilayer<=100; ilayer++)
    {
      TString layer("layer"); layer+= ilayer;
      for(int isd=0; isd<=2; isd++)
	{
	  TString sd("sd"); sd+= isd;
	  TString sdName("EE");
	  if(isd==1) sdName="HEfront";
	  if(isd==2) sdName="HEback";
	  
	  for(int idist=0; idist<sizeof(dists)/sizeof(TString); idist++)
	    {
	      TString dist=dists[idist];
	      c->Clear();
	      TH2F *h=(TH2F *)_file0->Get("analysis/"+dist+"_"+layer+"_"+sd);
	      if(h==0) continue;
	      h->Draw("colz");
	      if(dist.Contains("cell")) h->GetZaxis()->SetRangeUser(0,h->GetMaximum());
	      _file0->Get("analysis/boundary_"+layer+"_"+sd)->Draw("same");

	      h->GetYaxis()->SetTitleOffset(1.4);
	      h->GetZaxis()->SetTitleOffset(-0.5);

	      TPaveText *pt=new TPaveText(0.12,0.95,0.8,0.99,"brNDC");
	      pt->SetBorderSize(1);
	      pt->SetFillStyle(0);
	      pt->SetTextFont(42);
	      pt->SetTextAlign(12);
	      TString title(sdName + ", layer #"); title += ilayer; 
	      pt->AddText(title);
	      pt->Draw();

	      c->SaveAs(outDir+"/testhgcgeo_"+dist+"_"+layer+"_"+sd+".png");
	    }
	}
    }
}
