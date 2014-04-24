#!/usr/bin/env python

import ROOT
import sys
from ROOT import *
from UserCode.HGCanalysis.PlotUtils import *

def testSimHitsPos(url):

  customROOTstyle()
  gROOT.SetBatch(False)
  gStyle.SetPalette(51)

  fin=TFile.Open(url)
  HGC=fin.Get('hgcSimHitsAnalyzer/HGC')

  #auxiliary functions
  hit_phi='TMath::ATan2(hit_gy,hit_gx)'
  hit_rho='TMath::Sqrt(hit_gx*hit_gx+hit_gy*hit_gy+hit_gz*hit_gz)'
  hit_eta='0.5*TMath::Log((%s+hit_gz)/(%s-hit_gz))'%(hit_rho,hit_rho)
  deltaPhi='TVector2::Phi_mpi_pi({0}-{1})'
  deltaEta='{0}-{1}'

  #baseline cut
  cut="hit_type=={0}"

  #variables to plot
  vars={
    'dphi'    :[deltaPhi.format(hit_phi,'gen_phi[0]'), '#phi(hit)-#phi(gen) [rad];#hits', cut,                 ''     ,False, True,  False],
    'deta'    :[deltaEta.format(hit_eta,'gen_eta[0]'), '#eta(hit)-#eta(gen);#hits',       cut,                 ''     ,False, True,  False],
    'en'      :['hit_edep*1e3',                        'Energy [MeV];#hits',              cut,                 ''     ,False, True,  False],
    'xyzen'   :['hit_gy:hit_gx:hit_gz:hit_edep*1000',  'z [mm];x [mm];y [mm]; E [MeV]',   cut,                 'z'    ,False, False, False],
    'xz'      :['hit_gz:hit_gx',                       'x [mm];z[mm]',                    cut,                 ''     ,False, False, False],
    'yz'      :['hit_gz:hit_gy',                       'y [mm];z[mm]',                    cut,                 ''     ,False, False, False],
    'etaphi'  :[hit_eta+':'+hit_phi,                   '#phi [rad];#eta',                 cut,                 ''     ,False, False, False],
    'etaphien':[hit_eta+':'+hit_phi,                   '#phi [rad];#eta;E [MeV]',         'hit_edep*(%s)'%cut, 'colz' ,False, False, True]
    }     

  allCanvas=[]
  ihisto=0
  for v in vars:
    
    theVar=vars[v][0]
    theTitle=vars[v][1].split(';')
    theDrawOpt=vars[v][3]

    ic=len(allCanvas)
    allCanvas.append( ROOT.TCanvas("c"+v,"c"+v,1200,400) )
    allCanvas[ic].Divide(3,1);

    for i in xrange(0,3):

      theCut=vars[v][2].format(i)
      
      pad=allCanvas[ic].cd(i+1)
      pad.SetLogx(vars[v][4])
      pad.SetLogy(vars[v][5])
      pad.SetLogz(vars[v][6])
      if theDrawOpt.find('z')>0 : pad.SetRightMargin(0.2)
      HGC.Draw(theVar+'>>htemp',theCut,theDrawOpt)
      h=pad.GetPrimitive('htemp')
      try:
        ihisto=ihisto+1
        h.SetName('histo%d'%ihisto)
        if len(theTitle)>0 : h.GetXaxis().SetTitle(theTitle[0])
        if len(theTitle)>1 : h.GetYaxis().SetTitle(theTitle[1])
        if len(theTitle)>2 : h.GetZaxis().SetTitle(theTitle[2])
        h.GetXaxis().SetTitleSize(0.05)
        h.GetXaxis().SetLabelSize(0.04)
        h.GetXaxis().SetTitleOffset(1.2)
        h.GetYaxis().SetTitleSize(0.05)
        h.GetYaxis().SetLabelSize(0.04)
        h.GetYaxis().SetTitleOffset(1.2)
        h.GetZaxis().SetTitleSize(0.05)
        h.GetZaxis().SetLabelSize(0.04)
        h.GetZaxis().SetTitleOffset(1.2)
      except:
        print 'Failed for %s'%theVar
        
      if i==0: MyPaveText('CMS simulation')
      MyPaveText('[hit_type=%d]'%i,0.8,0.95,0.9,0.99).SetTextFont(42)

    allCanvas[ic].Modified()
    allCanvas[ic].Update()
    allCanvas[ic].SaveAs(v+'.png')


"""
checks the input arguments and steers the analysis
"""
def main():
  testSimHitsPos(sys.argv[1])

if __name__ == "__main__":
    main()
