{
  TH1F *pt=new TH1F("pth","pt",500,0,500);
  Float_t ptBins[]={5,10,15,20,30,50,75,100};
  for(size_t i=0; i<sizeof(ptBins)/sizeof(Float_t); i++) pt->Fill(ptBins[i]);
  TGraph *ptgr=new TGraph(pt);
  ptgr->SetName("pt");

  TH1F *rapidity=new TH1F("rap","rap",100,0,4);
  rapidity->SetName("rapidityh");
  Float_t yBins[]={1.5,1.75,2.0,2.25,2.5,2.75,2.9};
  for(size_t i=0; i<sizeof(yBins)/sizeof(Float_t); i++) rapidity->Fill(yBins[i]);
  TGraph *ygr=new TGraph(rapidity);
  ygr->SetName("rapidity");
  
  TFile *fOut=TFile::Open("particle_gun_pdf.root","RECREATE");
  ptgr->Write();
  ygr->Write();
  fOut->Close();
}
