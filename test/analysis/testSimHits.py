#!/usr/bin/env python

import ROOT
import sys
import optparse
import commands
from ROOT import *

from UserCode.HGCanalysis.PlotUtils import *
from UserCode.HGCanalysis.PileupUtils import *
from UserCode.HGCanalysis.HGCLayerUtils import *

"""
produces simple validation plots of the position of the sim hits for a particle gun sample
"""
def testSimHitsPos(fInUrl,accMap,sd,treeName='hgcSimHitsAnalyzer/HGC'):

  customROOTstyle()
  gROOT.SetBatch(False)
  gStyle.SetPalette(51)

  fin=TFile.Open(fInUrl)
  HGC=fin.Get(treeName)

  #loop over events
  for iev in xrange(0,HGC.GetEntriesFast()):
    HGC.GetEntry(iev)

    #require 1 single particle
    if HGC.ngen!=1 : continue
    genEta=HGC.gen_eta[0]
    genPhi=HGC.gen_phi[0]

    #require some hits
    if HGC.nhits==0 : continue
    for n in xrange(0,HGC.nhits):
      sdType=HGC.hit_type[n]
      if sdType!=sd : continue
      
      layer=HGC.hit_layer[n]
      sec=HGC.hit_sec[n]
      bin=HGC.hit_bin[n]
      rho,eta,phi=accMap[layer].getGlobalCoordinates(sec,bin)

      deltaPhi=TVector2.Phi_mpi_pi(phi-genPhi)
      deltaEta=eta-genEta

      print genPhi,deltaPhi,genEta,deltaEta

  fin.Close()

#  allCanvas=[]
#  ihisto=0
#  for v in vars:
#    
#    theVar=vars[v][0]
#    theTitle=vars[v][1].split(';')
#    theDrawOpt=vars[v][3]
#
#    ic=len(allCanvas)
#    allCanvas.append( ROOT.TCanvas("c"+v,"c"+v,1200,400) )
#    allCanvas[ic].Divide(3,1);
#
#    for i in xrange(0,3):
#
#      theCut=vars[v][2].format(i)
#      
#      pad=allCanvas[ic].cd(i+1)
#      pad.SetLogx(vars[v][4])
#      pad.SetLogy(vars[v][5])
#      pad.SetLogz(vars[v][6])
#      if theDrawOpt.find('z')>0 : pad.SetRightMargin(0.2)
#      HGC.Draw(theVar+'>>htemp',theCut,theDrawOpt)
#      h=pad.GetPrimitive('htemp')
#      try:
#        ihisto=ihisto+1
#        h.SetName('histo%d'%ihisto)
#        if len(theTitle)>0 : h.GetXaxis().SetTitle(theTitle[0])
#        if len(theTitle)>1 : h.GetYaxis().SetTitle(theTitle[1])
#        if len(theTitle)>2 : h.GetZaxis().SetTitle(theTitle[2])
#        h.GetXaxis().SetTitleSize(0.05)
#        h.GetXaxis().SetLabelSize(0.04)
#        h.GetXaxis().SetTitleOffset(1.2)
#        h.GetYaxis().SetTitleSize(0.05)
#        h.GetYaxis().SetLabelSize(0.04)
#        h.GetYaxis().SetTitleOffset(1.2)
#        h.GetZaxis().SetTitleSize(0.05)
#        h.GetZaxis().SetLabelSize(0.04)
#        h.GetZaxis().SetTitleOffset(1.2)
#      except:
#        print 'Failed for %s'%theVar
#        
#      if i==0: MyPaveText('CMS simulation')
#      MyPaveText('[hit_type=%d]'%i,0.8,0.95,0.9,0.99).SetTextFont(42)
#
#    allCanvas[ic].Modified()
#    allCanvas[ic].Update()
#    allCanvas[ic].SaveAs(v+'.png')


"""
checks the input arguments and steers the analysis
"""
def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in'         ,    dest='input'              , help='Input file',                          default=None           )
    parser.add_option('-s', '--sd'         ,    dest='sd'                 , help='Sensitive detector to analyse',       default=0,    type=int )
    (opt, args) = parser.parse_args()

    if opt.input is None :
        parser.print_help()
        sys.exit(1)
        
    #plot formatting
    customROOTstyle()
    #gROOT.SetBatch(True)
    #gROOT.SetBatch(False)
    gStyle.SetPalette(55)
    
    accMap=readSectorHistogramsFrom(fInUrl=opt.input,sd=opt.sd)
    testSimHitsPos(fInUrl=opt.input,accMap=accMap,sd=opt.sd)

if __name__ == "__main__":
    main()
