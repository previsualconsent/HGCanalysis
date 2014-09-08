{
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TString dists[]={"simcell","reccell","dx","dy","recix","reciy","simix","simiy"};
  TString layers[]={"layer1","layer2","layer3","layer4","layer5","layer6","layer7","layer8","layer9"};

  TFile *_file0 = TFile::Open("/tmp/psilva/HGCGeometry.root");

  TCanvas *c=new TCanvas("c","c",1500,500);

  for(int idist=0; idist<sizeof(dists)/sizeof(TString); idist++)
    {
      TString dist=dists[idist];
      for(int ilay=0; ilay<sizeof(layers)/sizeof(TString); ilay++)
	{
	  TString layer=layers[ilay];
	  
	  c->Clear();
	  c->Divide(3,1);
	  c->cd(1);
	  TH2F *h=(TH2F *)_file0->Get("analysis/"+dist+"_"+layer+"_sd0");
	  if(h){
	    h->Draw("colz");
	    if(dist.Contains("cell")) h->GetZaxis()->SetRangeUser(0,h->GetMaximum());
	    _file0->Get("analysis/boundary_"+layer+"_sd0")->Draw("same");
	  }
	  c->cd(2);
	  TH2F *h=(TH2F *)_file0->Get("analysis/"+dist+"_"+layer+"_sd1");
	  if(h){
	    h->Draw("colz");
	    if(dist.Contains("cell")) h->GetZaxis()->SetRangeUser(0,h->GetMaximum());
	    _file0->Get("analysis/boundary_"+layer+"_sd1")->Draw("same");
	  }
	  c->cd(3);
	  TH2F *h=(TH2F *)_file0->Get("analysis/"+dist+"_"+layer+"_sd2");
	  if(h){
	    h->Draw("colz");
	    if(dist.Contains("cell")) h->GetZaxis()->SetRangeUser(0,h->GetMaximum());
	    _file0->Get("analysis/boundary_"+layer+"_sd2")->Draw("same");
	  }

	  c->SaveAs("~/www/testhgcgeo_"+dist+"_"+layer+".png");
	  //c->SaveAs("~/www/testhgcgeo_"+dist+"_"+layer+".C");
	}
    }
}
