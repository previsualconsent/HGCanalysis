{
  TH1F *pt=new TH1F("pth","pt",500,0,500);
  Float_t ptBins[]={2.5,5,10,15,20,30,50,60,80};
  for(size_t i=0; i<sizeof(ptBins)/sizeof(Float_t); i++) pt->Fill(ptBins[i]);
  TGraph *ptgr=new TGraph(pt);
  ptgr->SetName("pt");

  TH1F *rapidity=new TH1F("rap","rap",100,0,4);
  rapidity->SetName("rapidityh");
  Float_t yBins[]={1.5,1.7,1.9,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9};
  for(size_t i=0; i<sizeof(yBins)/sizeof(Float_t); i++) rapidity->Fill(yBins[i]);
  TGraph *ygr=new TGraph(rapidity);
  ygr->SetName("rapidity");
  
  TFile *fOut=TFile::Open("particle_gun_pdf.root","RECREATE");
  ptgr->Write();
  ygr->Write();
  fOut->Close();
}
