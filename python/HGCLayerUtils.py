#!/usr/bin/env python

import ROOT
from ROOT import *

"""
Wrapper for an HGCal
"""
class LayerAccumulator:
    
    def __init__(self,baseHisto,layer):

        #clone and rebin to new cell width
        self.layer=layer
        self.accumulator={}
        self.accumulator[0]=baseHisto.Clone('E_accumulator_%d_0'%layer)
        self.accumulator[0].SetDirectory(0)
        self.accumulator[0].Reset('ICE')
                        
        #create the eta map
        newName=self.accumulator[0].GetName()
        self.rhoMap=self.accumulator[0].Clone(newName.replace('E_','rho_'))
        self.rhoMap.SetDirectory(0)
        self.etaMap=self.accumulator[0].Clone(newName.replace('E_','eta_'))
        self.etaMap.SetDirectory(0)
        self.phiMap=self.accumulator[0].Clone(newName.replace('E_','phi_'))
        self.phiMap.SetDirectory(0)

    def defineCoordinates(self, gxH, gyH, gzH):
        
        #fill the eta map with the coordinates
        for xbin in xrange(1,self.accumulator[0].GetXaxis().GetNbins()+1):
            for ybin in xrange(1,self.accumulator[0].GetYaxis().GetNbins()+1):
                x=gxH.GetBinContent(xbin,ybin)
                y=gyH.GetBinContent(xbin,ybin)
                z=gzH.GetBinContent(xbin,ybin)
                rho=TMath.Sqrt(x*x+y*y+z*z)
                phi=TMath.ATan2(y,x)
                eta=0
                if rho>z: eta=0.5*TMath.Log( (rho+z)/(rho-z) )
                self.rhoMap.SetBinContent(xbin,ybin,rho)
                self.etaMap.SetBinContent(xbin,ybin,eta)
                self.phiMap.SetBinContent(xbin,ybin,phi)
                
    def getGlobalCoordinates(self,sector,ibin):
        xbin,ybin,zbin = ROOT.Long(), ROOT.Long(), ROOT.Long()
        self.rhoMap.GetBinXYZ(ibin,xbin,ybin,zbin)
        return self.rhoMap.GetBinContent(xbin,ybin),self.etaMap.GetBinContent(xbin,ybin),self.phiMap.GetBinContent(xbin,ybin)
                    
    def addSector(self, sector):
        if sector==0: return
        self.accumulator[sector]=self.accumulator[0].Clone('E_accumulator_%d_%d'%(self.layer,sector)) 
        self.accumulator[sector].SetDirectory(0)

    def reset(self) :
        for sector in self.accumulator:
            self.accumulator[sector].Reset('ICE')

    def fill(self,edep,ibin,sector):
        self.accumulator[sector].SetBinContent(ibin,edep+self.accumulator[sector].GetBinContent(ibin))

                            
"""
get all sector histograms into memory
"""
def readSectorHistogramsFrom(fInUrl,sd,baseDir='hgcSimHitsAnalyzer'):
    accumulatorsMap={}

    #read energy accumulators from map and define the coordinates
    fIn=ROOT.TFile.Open(fInUrl)
    dir=fIn.Get(baseDir)
    for key in dir.GetListOfKeys():
        if key.GetName().find('E_')<0 : continue
        keyName=key.GetName().replace('E_','')
        sectorInfo=keyName.split('_')
        if int(sectorInfo[1])!=sd : continue
        layer=int(sectorInfo[2])
        sector=int(sectorInfo[3])
        if layer in accumulatorsMap:
            accumulatorsMap[layer].addSector(sector)
        else:
            accumulatorsMap[ layer ] = LayerAccumulator( fIn.Get(baseDir+'/E_'+keyName), layer )
            accumulatorsMap[ layer ].defineCoordinates( fIn.Get(baseDir+'/gx_'+keyName), fIn.Get(baseDir+'/gy_'+keyName), fIn.Get(baseDir+'/gz_'+keyName) )
    fIn.Close()

    #all done here
    for layer in accumulatorsMap:
        print 'Layer=%d has %d sectors'%(layer,len(accumulatorsMap[layer].accumulator))
    return accumulatorsMap

