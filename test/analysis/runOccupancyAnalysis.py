#!/usr/bin/env python

import ROOT
import sys
import optparse
import commands
from array import array
import numpy as np
from ROOT import *
from UserCode.HGCanalysis.PlotUtils import *
from PileupUtils import *
                            

"""
define a subsector histogram for a given layer based on its dimensions and cell size
"""
def getSubsectorTemplateForLayer(layer,h,bl,tl,cell):
    ndivx=int(2*tl/cell)
    ndivy=int(2*h/cell)
    name='layer%dcell%d'%(layer,cell)
    title='Layer %d Cell %dx%d mm^{2};x [mm];y [mm]'%(layer,cell,cell)
    h=TH2F(name,title,ndivx,-tl-cell/2,tl-cell/2,ndivy,-h-cell/2,h-cell/2)
    h.SetDirectory(0)
    return h

"""
loops of the layers and instantiates the subsector histograms
"""
def getLayerSubSectorHistos(fin,cell):

    #get geometry definitions from file
    eeHalfHeight      = fin.Get('hgcSimHitsAnalyzer/eeHalfHeight')
    eeBottomHalfWidth = fin.Get('hgcSimHitsAnalyzer/eeBottomHalfWidth')
    eeTopHalfWidth    = fin.Get('hgcSimHitsAnalyzer/eeTopHalfWidth')

    #instantiate histograms
    subsectorHistos={}
    for xbin in xrange(1,eeHalfHeight.GetXaxis().GetNbins()):
        layer=int(eeHalfHeight.GetXaxis().GetBinLowEdge(xbin))
        hl=eeHalfHeight.GetBinContent(xbin)
        bl=eeBottomHalfWidth.GetBinContent(xbin)
        tl=eeTopHalfWidth.GetBinContent(xbin)
        if hl==0 : continue
        templH=getSubsectorTemplateForLayer(layer,hl,bl,tl,cell)
        subsectorHistos[layer]=[]
        for isec in xrange(0,18):
            h=templH.Clone('%s_sec%d'%(templH.GetName(),isec+1))
            h.SetDirectory(0)
            subsectorHistos[layer].append( h )

    #all done here
    return subsectorHistos


"""
loops over the events and collects the electron energies
"""
def runOccupancyAnalysis(url='HGCSimHitsAnalysis.root',avgPU=0,mipEn=54.8,cellSize=10,treeName='hgcSimHitsAnalyzer/HGC') :

    customROOTstyle()
    gROOT.SetBatch(True)
    #gROOT.SetBatch(False)
    gStyle.SetPalette(55)

    #analyze generated events
    fin=TFile.Open(url)
    Events=fin.Get(treeName)

    #base templates in global coordinates and in local coordinates
    layerHistos=getLayerSubSectorHistos(fin,cellSize)
    etaArray=array('d', sorted([-5]+np.arange(-3.4, -1.2, 0.05,dtype=np.double).tolist()+[-1.5,0,1.5]+np.arange(1.2,3.4, 0.05,dtype=np.double).tolist()+[5]) )
    etaPhiTemplate=TH2F('etaphi',';Pseudo-rapidity;#phi [rad];Energy',len(etaArray)-1,etaArray,64,-3.2,3.2)

    #deposits in sectors are mapped in local coordinates
    layerEnAverage={}
    layerEnAverageEtaPhi={}
    layerOccAverage1mip={}
    layerOccAverage10mip={}
    layerOccAverage20mip={}
    layerYtoEtaMap={}                                   #*central* y to eta translation
    eeTranslZ=fin.Get('hgcSimHitsAnalyzer/eeTranslZ')
    eeTranslY=fin.Get('hgcSimHitsAnalyzer/eeTranslY')

    #create histograms
    for l in layerHistos:

        #for energy averages
        layerEnAverage[l]=layerHistos[l][0].Clone('avgenergy_layer%d'%l)
        layerEnAverage[l].SetDirectory(0)
        layerEnAverageEtaPhi[l]=etaPhiTemplate.Clone('avgenergy_etaphi_layer%d'%l)
        layerEnAverageEtaPhi[l].SetDirectory(0)

        #for ocuppancy averages
        layerOccAverage1mip[l]=layerHistos[l][0].Clone('avgocc_1mip_layer%d'%l)
        layerOccAverage1mip[l].SetDirectory(0)
        layerOccAverage10mip[l]=layerHistos[l][0].Clone('avgocc_10mip_layer%d'%l)
        layerOccAverage10mip[l].SetDirectory(0)
        layerOccAverage20mip[l]=layerHistos[l][0].Clone('avgocc_20mip_layer%d'%l)
        layerOccAverage20mip[l].SetDirectory(0)

        #to map y to eta
        ytoEtaMap=TH1F('ytoetamap_layer%d'%l,';y @ x=0 [mm]; Pseudo-rapidity',
                       layerHistos[l][0].GetYaxis().GetNbins(),
                       layerHistos[l][0].GetYaxis().GetXmin(),
                       layerHistos[l][0].GetYaxis().GetXmax())
        ytoEtaMap.SetDirectory(0)
        for xbin in xrange(1,layerHistos[l][0].GetYaxis().GetNbins()):
            y=layerHistos[l][0].GetYaxis().GetBinCenter(xbin)+eeTranslY.GetBinContent(l)
            z=eeTranslZ.GetBinContent(l)
            p4=TLorentzVector(0,y,z,TMath.Sqrt(y*y+z*z))
            ytoEtaMap.SetBinContent(xbin,p4.Eta())
        layerYtoEtaMap[l]=ytoEtaMap

        
    #generate pileup sets
    evRanges=generateInTimePUSets(Events.GetEntriesFast(),avgPU)
    
    #iterate over events
    iSetCtr=0
    for iSet in evRanges:

        iSetCtr=iSetCtr+1
                
        #reset entries in current histos
        for l in layerHistos:
            for h in layerHistos[l]:
                h.Reset('ICE')
        
        #sum up N events
        for iev in xrange(iSet[0],iSet[1]+1):

            #get hits for new event
            Events.GetEntry(iev)
            if Events.nee==0 : continue

            #get energy deposits info for ECAL
            for idep in xrange(0,Events.nee):

                edep=Events.ee_edep[idep]*1e6/mipEn
                layer=Events.ee_layer[idep]

                #for average energy deposits per layer in *global* coordinates
                gx=Events.ee_gx[idep]
                gy=Events.ee_gy[idep]
                gz=Events.ee_gz[idep]
                radius=TMath.Sqrt(gx*gx+gy*gy+gz*gz)
                p4=TLorentzVector(gx,gy,gz,radius)
                phi=p4.Phi()
                eta=p4.Eta()
                layerEnAverageEtaPhi[layer].Fill(eta,phi,edep)
                
                #only positive 
                if Events.ee_zp[idep]<1: continue

                #fill sector histos
                x=Events.ee_x[idep]*Events.ee_subsec[idep]
                y=Events.ee_y[idep]
                layerHistos[ layer ][ Events.ee_sec[idep]-1 ].Fill(x,y,edep)
                layerEnAverage[ layer ].Fill(x,y,edep/len(layerHistos[layer]))
                
        #now analyze
        for l in layerHistos:

            #occupancy
            for secH in layerHistos[l]:
                for xbin in xrange(1,secH.GetXaxis().GetNbins()+1) :
                    for ybin in xrange(1,secH.GetYaxis().GetNbins()+1) :
                        enInCell=secH.GetBinContent(xbin,ybin)
                        if enInCell<1: continue
                        curOcc=layerOccAverage1mip[l].GetBinContent(xbin,ybin)+1./len(layerHistos[l])
                        layerOccAverage1mip[l].SetBinContent(xbin,ybin,curOcc)
                        if enInCell<10: continue
                        curOcc=layerOccAverage10mip[l].GetBinContent(xbin,ybin)+1./len(layerHistos[l])
                        layerOccAverage10mip[l].SetBinContent(xbin,ybin,curOcc)
                        if enInCell<20: continue
                        curOcc=layerOccAverage20mip[l].GetBinContent(xbin,ybin)+1./len(layerHistos[l])
                        layerOccAverage20mip[l].SetBinContent(xbin,ybin,curOcc)
                    
    #save profiles obtained as function of y and y->eta map
    fout=TFile.Open('occprofiles_%d.root'%avgPU,'RECREATE')
    for l in layerHistos:
        layerEnAverage[l].SetDirectory(fout)
        layerEnAverage[l].Scale(1./iSetCtr)
        layerEnAverage[l].Write()
        layerEnAverageEtaPhi[l].SetDirectory(fout)
        layerEnAverageEtaPhi[l].Scale(1./iSetCtr)
        layerEnAverageEtaPhi[l].Write()
        layerYtoEtaMap[l].SetDirectory(fout)
        layerYtoEtaMap[l].Write()
        layerOccAverage1mip[l].Scale(1./iSetCtr)
        layerOccAverage1mip[l].SetDirectory(fout)
        layerOccAverage1mip[l].Write()
        layerOccAverage10mip[l].Scale(1./iSetCtr)
        layerOccAverage10mip[l].SetDirectory(fout)
        layerOccAverage10mip[l].Write()
        layerOccAverage20mip[l].Scale(1./iSetCtr)
        layerOccAverage20mip[l].SetDirectory(fout)
        layerOccAverage20mip[l].Write()
    fout.Close()

    
    fin.Close()
    



"""
checks the input arguments and steers the analysis
"""
def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in'         ,    dest='input'              , help='Input file', default=None)
    parser.add_option('-p', '--pu'         ,    dest='pu'                 , help='Average PU to overlay', default=0, type=int)
    (opt, args) = parser.parse_args()

    if opt.input is None :
        parser.print_help()
        sys.exit(1)

    print 'Overlaying %d average MinBias events'%opt.pu

    #later can add other options
    runOccupancyAnalysis(opt.input,opt.pu)

if __name__ == "__main__":
    main()


