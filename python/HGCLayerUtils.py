#!/usr/bin/env python

import ROOT
from ROOT import *

"""
Useful to integrate in the vicinity of a cell with different granularities
"""
class CellIntegrator:
    def __init__(self):
        self.center      = {1:[2,2],  2:[4,4],  3:[5,5]}
        self.step        = {1:3,      2:6,      3:9}
        self.integRange  = {1:[-1,1], 2:[-3,2], 3:[-4,4]}
        self.signalRange = {1:[0,0],  2:[-1,0], 3:[-1,1]}
    def __isSignal(self,xbin,ybin,cell):
        if xbin<self.signalRange[cell][0] or xbin>self.signalRange[cell][1]: return False
        if ybin<self.signalRange[cell][0] or ybin>self.signalRange[cell][1]: return False
        return True
    def getSRNumber(self,xbin,ybin,cell):
        if self.__isSignal(xbin,ybin,cell): return 0
        elif cell==1:
            if xbin<0 or ybin<0           : return 2
            else                          : return 1
        elif cell==2:
            if xbin<-1 or ybin<-1         : return 2
            else                          : return 1
        elif cell==3:
            if xbin<-1 or ybin<-1         : return 2
            else                          : return 1

        

"""
Wrapper for an HGCal layer
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
        self.rhoMap={}
        self.rhoMap[0]=self.accumulator[0].Clone(newName.replace('E_','rho_'))
        self.rhoMap[0].SetDirectory(0)
        self.etaMap={}
        self.etaMap[0]=self.accumulator[0].Clone(newName.replace('E_','eta_'))
        self.etaMap[0].SetDirectory(0)
        self.phiMap={}
        self.phiMap[0]=self.accumulator[0].Clone(newName.replace('E_','phi_'))
        self.phiMap[0].SetDirectory(0)

    def defineCoordinates(self, sector, gxH, gyH, gzH):
        
        #fill the eta map with the coordinates
        for xbin in xrange(1,self.accumulator[sector].GetXaxis().GetNbins()+1):
            for ybin in xrange(1,self.accumulator[sector].GetYaxis().GetNbins()+1):
                x=gxH.GetBinContent(xbin,ybin)
                y=gyH.GetBinContent(xbin,ybin)
                z=gzH.GetBinContent(xbin,ybin)
                rho=TMath.Sqrt(x*x+y*y+z*z)

                #skip masked channels (0 netry)
                if rho==0 or z==0 : continue

                phi=TMath.ATan2(y,x)
                eta=0
                if rho>z: eta=0.5*TMath.Log( (rho+z)/(rho-z) )
                self.rhoMap[sector].SetBinContent(xbin,ybin,rho)
                self.etaMap[sector].SetBinContent(xbin,ybin,eta)
                self.phiMap[sector].SetBinContent(xbin,ybin,phi)
        
    def getGlobalCoordinates(self,sector,ibin):
        xbin,ybin,zbin = ROOT.Long(), ROOT.Long(), ROOT.Long()
        self.rhoMap[sector].GetBinXYZ(ibin,xbin,ybin,zbin)
        return self.rhoMap[sector].GetBinContent(xbin,ybin),self.etaMap[sector].GetBinContent(xbin,ybin),self.phiMap[sector].GetBinContent(xbin,ybin)
                    
    def addSector(self, sector, gxH, gyH, gzH):
        #clone first sector if needed
        if not(sector in self.accumulator):
            self.accumulator[sector]=self.accumulator[0].Clone('E_accumulator_%d_%d'%(self.layer,sector)) 
            self.accumulator[sector].SetDirectory(0)
            self.rhoMap[sector]=self.rhoMap[0].Clone('rho_accumulator_%d_%d'%(self.layer,sector))
            self.rhoMap[sector].SetDirectory(0)
            self.rhoMap[sector].Reset('ICE')
            self.etaMap[sector]=self.etaMap[0].Clone('eta_accumulator_%d_%d'%(self.layer,sector))
            self.etaMap[sector].SetDirectory(0)
            self.etaMap[sector].Reset('ICE')
            self.phiMap[sector]=self.phiMap[0].Clone('phi_accumulator_%d_%d'%(self.layer,sector))
            self.phiMap[sector].SetDirectory(0)
            self.phiMap[sector].Reset('ICE')
        self.defineCoordinates(sector,gxH,gyH,gzH)

    def reset(self) :
        for sector in self.accumulator:
            self.accumulator[sector].Reset('ICE')

    def fill(self,edep,ibin,sector):
        try:
            return self.accumulator[sector].SetBinContent(ibin,edep+self.accumulator[sector].GetBinContent(ibin))
        except:
            return 0
                            
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
        if not(layer in accumulatorsMap):
            accumulatorsMap[ layer ] = LayerAccumulator( fIn.Get(baseDir+'/E_'+keyName), layer )
        accumulatorsMap[layer].addSector( sector, fIn.Get(baseDir+'/gx_'+keyName) ,fIn.Get(baseDir+'/gy_'+keyName), fIn.Get(baseDir+'/gz_'+keyName) )
        
    fIn.Close()

    #all done here
    for layer in accumulatorsMap:
        print 'Layer=%d has %d sectors'%(layer,len(accumulatorsMap[layer].accumulator))
    return accumulatorsMap

