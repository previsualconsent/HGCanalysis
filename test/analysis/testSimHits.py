#!/usr/bin/env python

import ROOT
import sys
import optparse
import commands
import array
from ROOT import *

from UserCode.HGCanalysis.PlotUtils import *
from UserCode.HGCanalysis.PileupUtils import *
from UserCode.HGCanalysis.HGCLayerUtils import *

"""
Wrapper to store hits
"""
class HitIntegrator:
  def __init__(self,accMap):
    self.subDets=[]
    self.varNames=[]
    for subDet in accMap:
      for layer in accMap[subDet]:
        self.subDets.append( self.getKeyFor(subDet, layer ) )
    self.vars=['avgdeta','avgdphi','en','endr01']
    self.hitsCollection={}
    for subDet in self.subDets:
      self.hitsCollection[subDet]={}
      for var in self.vars:
        self.hitsCollection[subDet][var]=0
        self.varNames.append('%s_%s'%(subDet,var))

  def getKeyFor(self,subDet,layer):
    key='s%d'%subDet
    if layer<0 : key+='m'
    else       : key+='p'
    key+='%d'%abs(layer)
    return key
  
  def resetHits(self):
    for subDet in self.hitsCollection:
      for var in self.hitsCollection[subDet]:
        self.hitsCollection[subDet][var]=0

  def getHits(self)     : return self.hitsCollection
  def getVarNames(self) : return self.varNames
  def getVarVals(self)  :
    varVals=[]
    for ivar in self.varNames:
      subDet,var=ivar.split('_')
      varVals.append( self.hitsCollection[subDet][var] )
    return varVals
  
  def integrateHit(self,subDet,layer,edep,deltaEta,deltaPhi):

    #assign sub detector
    key=self.getKeyFor(subDet,layer)

    #compute the variables of interest
    deltaR=TMath.Sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi)
    iVarVal={'avgdeta'    : deltaEta*edep,
             'avgdphi'    : deltaPhi*edep,
             'en'       : edep,
             'endr01'  : 0. if deltaR>0.1  else edep
             }
    for var in iVarVal:
      curVal=self.hitsCollection[key][var]+iVarVal[var]
      self.hitsCollection[key][var]=curVal
    

"""
produces simple validation plots of the position of the sim hits for a particle gun sample
"""
def integrateSimHits(fInUrl,fout,accMap,treeName='hgcSimHitsAnalyzer/HGC'):

  customROOTstyle()
  gROOT.SetBatch(False)
  gStyle.SetPalette(51)

  hitIntegrator=HitIntegrator(accMap)

  fout.cd()
  output_tuple = ROOT.TNtuple("HGC","HGC","genPt:genEta:genPhi:"+":".join(hitIntegrator.getVarNames()))
  output_tuple.SetDirectory(fout)
  
  #loop over events
  HGC=ROOT.TChain(treeName)
  HGC.Add(fInUrl)
  for iev in xrange(0,HGC.GetEntries()):
    HGC.GetEntry(iev)

    sys.stdout.write( '\r Status [%d/%d]'%(iev,HGC.GetEntries()))
       
    if HGC.ngen==0 or HGC.ngen>2 : continue
    if HGC.nhits==0 : continue

    #create a new entry for each new particle (particle gun may have two shot in opposite eta)
    for ngen in xrange(0,HGC.ngen):

      genPt=HGC.gen_pt[ngen]
      genEta=HGC.gen_eta[ngen]
      genPhi=HGC.gen_phi[ngen]

      hitIntegrator.resetHits()

      for n in xrange(0,HGC.nhits):
      
        sdType=HGC.hit_type[n]
        layer=HGC.hit_layer[n]
        sec=HGC.hit_sec[n]
        bin=HGC.hit_bin[n]
        rho,eta,phi=accMap[sdType][layer].getGlobalCoordinates(sec,bin)

        #filter hits to match the eta of the generated particle
        if eta*genEta<0 : continue

        #integrate hit
        edep=HGC.hit_edep[n]
        deltaPhi=TVector2.Phi_mpi_pi(phi-genPhi)
        deltaEta=eta-genEta
        hitIntegrator.integrateHit(sdType,abs(layer),edep,deltaEta,deltaPhi)

      #fill the ntuple
      varVals=[genPt,genEta,genPhi]+hitIntegrator.getVarVals()
      output_tuple.Fill(array.array("f",varVals))
      if iev%5000 : output_tuple.AutoSave('SaveSelf')

  #write the ntuple
  fout.cd()
  output_tuple.Write()

"""
checks the input arguments and steers the analysis
"""
def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in' ,  dest='input', help='Input directory',    default=None)
    parser.add_option('-t', '--tag',  dest='tag'  , help='Files matching tag', default=None)
    (opt, args) = parser.parse_args()

    #check inputs
    if opt.input is None or opt.tag is None:
        parser.print_help()
        sys.exit(1)

    #filter input files
    lsOutput = commands.getstatusoutput('cmsLs %s | grep root | awk \'{print $5}\''%(opt.input))[1].split()
    fInUrl=[]
    for f in lsOutput:
      if f.find(opt.tag)<0 : continue
      fInUrl.append( commands.getstatusoutput('cmsPfn '+f)[1] )
    if len(fInUrl)==0:
      print 'No files matching %s in %s have been found'%(opt.tag,opt.input)
      parser.print_help()
      sys.exit(1)

    #read geometry
    accMap={ 0 : readSectorHistogramsFrom(fInUrl=fInUrl[0],sd=0),
             1 : readSectorHistogramsFrom(fInUrl=fInUrl[0],sd=1),
             2 : readSectorHistogramsFrom(fInUrl=fInUrl[0],sd=2) }

    #run analysis
    ifile=0
    for f in fInUrl:
      ifile=ifile+1
      print ifile,f
      fout = ROOT.TFile("%s_%d.root"%(opt.tag,ifile), "RECREATE")
      integrateSimHits(fInUrl=f,fout=fout,accMap=accMap)
      fout.Close()

if __name__ == "__main__":
    main()
