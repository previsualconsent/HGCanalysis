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
    ndivx=int(bl/cell)
    ndivx=(ndivx+int((tl-ndivx*cell)/cell))
    ndivy=int(h/cell)
    name='layer%dcell%d'%(layer,cell)
    title='Layer %d Cell %dx%d mm^{2};x [mm];y [mm]'%(layer,cell,cell)
    h=TH2F(name,title,2*ndivx,-cell*ndivx,cell*ndivx,2*ndivy,-cell*ndivy,cell*ndivy)
    h.SetDirectory(0)
    return h

"""
loops of the layers and instantiates the subsector histograms
"""
def getLayerSubSectorHistos(fin,cell,sdType,sdCtrPeriod):

    #get geometry definitions from file
    sdHalfHeight      = fin.Get('hgcSimHitsAnalyzer/sens%d_HalfHeight'%sdType)
    sdBottomHalfWidth = fin.Get('hgcSimHitsAnalyzer/sens%d_BottomHalfWidth'%sdType)
    sdTopHalfWidth    = fin.Get('hgcSimHitsAnalyzer/sens%d_TopHalfWidth'%sdType)
    sdNSectors        = fin.Get('hgcSimHitsAnalyzer/sens%d_NSectors'%sdType)

    #instantiate histograms
    subsectorHistos={}
    sdCtr=0
    for xbin in xrange(1,sdHalfHeight.GetXaxis().GetNbins()):

        layer=int(sdHalfHeight.GetXaxis().GetBinLowEdge(xbin))
        hl=sdHalfHeight.GetBinContent(xbin)
        bl=sdBottomHalfWidth.GetBinContent(xbin)
        tl=sdTopHalfWidth.GetBinContent(xbin)

        if hl<1 : continue
        if layer<0: continue

        #create new histogram only for the given sample period
        if sdCtrPeriod==1 or sdCtr%sdCtrPeriod==0 :

            layer=layer/sdCtrPeriod+1

            templH=getSubsectorTemplateForLayer(layer,hl,bl,tl,cell)
            subsectorHistos[layer]=[]
    
            xbin=sdNSectors.GetXaxis().FindBin(layer)
            nsectors=int(sdNSectors.GetBinContent(xbin))
        
            for isec in xrange(0,nsectors):
                h=templH.Clone('%s_sec%d'%(templH.GetName(),isec+1))
                h.SetDirectory(0)
                subsectorHistos[layer].append( h )

        sdCtr=sdCtr+1

    #all done here
    return subsectorHistos

"""
loops over the events and collects the electron energies
"""
def runOccupancyAnalysis(url='HGCSimHitsAnalysis.root',avgPU=0,sdType=0,cellSize=10,mipEn=54.8,treeName='hgcSimHitsAnalyzer/HGC') :

    #analyze generated events
    fin=TFile.Open(url)
    Events=fin.Get(treeName)

    #base templates in global coordinates and in local coordinates
    sdCtrPeriod=1
    subLayerToSample=[0]
    if sdType in [0,1]:
        sdCtrPeriod=3
        subLayerToSample=[0,1]
    layerHistos=getLayerSubSectorHistos(fin,cellSize,sdType,sdCtrPeriod)

    #create histograms for the analysis
    layerOccAverage={}
    layerEnRanges={}
    layerAvgEnRanges={}
    for l in layerHistos:

        #occupancy profiles
        for thr in [1,5,10,25,50,100]:
            if not thr in layerOccAverage : layerOccAverage[thr]={}
            layerOccAverage[thr][l]=layerHistos[l][0].Clone('avgocc_%dmip_layer%d'%(thr,l))
            layerOccAverage[thr][l].SetDirectory(0)

        #energy ranges
        for q in ['max','mean','std','q50','q90']:
            if not q in layerEnRanges :
                layerEnRanges[q]={}
                layerAvgEnRanges[q]={}
            layerEnRanges[q][l]=layerHistos[l][0].Clone('eranges_%s_layer%d'%(q,l))
            layerEnRanges[q][l].SetDirectory(0)
            layerAvgEnRanges[q][l]=layerHistos[l][0].Clone('avgeranges_%s_layer%d'%(q,l))
            layerAvgEnRanges[q][l].SetDirectory(0)

    #dictionary to translate to pseudo-rapidity
    cellEtaMap={}
    sdTranslZ=fin.Get('hgcSimHitsAnalyzer/sens%d_TranslZ'%sdType)
    sdTranslX=fin.Get('hgcSimHitsAnalyzer/sens%d_TranslX'%sdType)
    for l in layerHistos:
        cellEtaMap[l]=TH1F('celletamap_layer%d'%l,';y @ x=0 [mm]; Pseudo-rapidity',
                           layerHistos[l][0].GetYaxis().GetNbins(),
                           layerHistos[l][0].GetYaxis().GetXmin(),
                           layerHistos[l][0].GetYaxis().GetXmax())
        cellEtaMap[l].SetDirectory(0)
        lbin=sdTranslX.GetXaxis().FindBin(l)
        for xbin in xrange(1,layerHistos[l][0].GetYaxis().GetNbins()):
            gx=layerHistos[l][0].GetYaxis().GetBinCenter(xbin)+sdTranslX.GetBinContent(lbin)
            gz=sdTranslZ.GetBinContent(lbin)
            g=TMath.Sqrt(gx*gx+gz*gz)
            eta=0.5*TMath.Log( (g+gz)/(g-gz) )
            cellEtaMap[l].SetBinContent(xbin,eta)
        
    #generate pileup sets
    evRanges=generateInTimePUSets(Events.GetEntriesFast(),avgPU)
    
    #iterate over events
    iSetCtr=0
    for iSet in evRanges:

        #if iSetCtr>1: break

        iSetCtr=iSetCtr+1
        sys.stdout.write( '\r Status [%d/%d]'%(iSetCtr,len(evRanges)) )
        sys.stdout.flush()
                
        #reset entries in current histos
        for l in layerHistos:
            for h in layerHistos[l]:
                h.Reset('ICE')
        
        #sum up N events
        for iev in xrange(iSet[0],iSet[1]+1):

            #get hits for new event
            Events.GetEntry(iev)
            if Events.nhits==0 : continue

            #get energy deposits info for the sensitive detector (only positive in z)
            for idep in xrange(0,Events.nhits):

                if Events.hit_type[idep]!=sdType : continue
                if Events.hit_zp[idep]<1         : continue
                
                edep=Events.hit_edep[idep]*1e6/mipEn
                layer=Events.hit_layer[idep]
                sector=Events.hit_sec[idep]-1
                x=Events.hit_x[idep]
                y=Events.hit_y[idep]

                subLayer=layer%sdCtrPeriod+1
                layer=layer/sdCtrPeriod
                if not subLayer in subLayerToSample: continue
                
                #fill base histos
                layerHistos[ layer ][ sector ].Fill(x,y,edep)

                
        #now analyze
        for l in layerHistos:
        
            for xbin in xrange(1,layerHistos[l][0].GetXaxis().GetNbins()+1) :
                for ybin in xrange(1,layerHistos[l][0].GetYaxis().GetNbins()+1) :

                    #occupancy
                    edeps=[]
                    for secH in layerHistos[l]:
                        
                        enInCell=secH.GetBinContent(xbin,ybin)
                        edeps.append(enInCell)
                        
                        for thr in layerOccAverage:
                            if enInCell<thr : continue
                            curOcc=layerOccAverage[thr][l].GetBinContent(xbin,ybin)+1./len(layerHistos[l])
                            layerOccAverage[thr][l].SetBinContent(xbin,ybin,curOcc)


                    #quantiles for the energy deposition
                    if max(edeps)==0 : continue
                    
                    for q in layerEnRanges:

                        val=0
                        if q=='max':  val=max(edeps)
                        if q=='mean': val=np.mean(edeps)
                        if q=='std':  val=np.std(edeps)
                        if q=='q90':  val=np.percentile(edeps,90)
                        if q=='q50':  val=np.percentile(edeps,50)
                        
                        maxVal=max(layerEnRanges[q][l].GetBinContent(xbin,ybin),val)
                        avgVal=layerAvgEnRanges[q][l].GetBinContent(xbin,ybin)+val/len(edeps)
                        
                        layerAvgEnRanges[q][l].SetBinContent(xbin,ybin,avgVal)
                        layerEnRanges[q][l].SetBinContent(xbin,ybin,maxVal)
                    
                    
                    
    #all done with the inputs
    fin.Close()

                    
    #save profiles obtained as function of y and y->eta map
    foutUrl='occprofiles_pu%d_sd%d_cell%d.root'%(avgPU,sdType,cellSize)
    fout=TFile.Open(foutUrl,'RECREATE')
    for l in layerHistos:
        cellEtaMap[l].SetDirectory(fout)
        cellEtaMap[l].Write()
        for thr in layerOccAverage:
            layerOccAverage[thr][l].Scale(1./iSetCtr)
            layerOccAverage[thr][l].SetDirectory(fout)
            layerOccAverage[thr][l].Write()
        for q in layerEnRanges:
            layerEnRanges[q][l].SetDirectory(fout)
            layerEnRanges[q][l].Write()
            layerAvgEnRanges[q][l].Scale(1./iSetCtr)
            layerAvgEnRanges[q][l].SetDirectory(fout)
            layerAvgEnRanges[q][l].Write()
                        
    fout.Close()
    print 'Results stored in %s'%foutUrl


"""
checks the input arguments and steers the analysis
"""
def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in'         ,    dest='input'              , help='Input file', default=None)
    parser.add_option('-p', '--pu'         ,    dest='pu'                 , help='Average PU to overlay', default=0, type=int)
    parser.add_option('-s', '--sd'         ,    dest='sd'                 , help='Sensitive detector to analyse', default=0, type=int)
    parser.add_option('-c', '--cell'       ,    dest='cell'               , help='cell size [mm]', default=10, type=int)
    (opt, args) = parser.parse_args()

    if opt.input is None :
        parser.print_help()
        sys.exit(1)


    #plot formatting
    customROOTstyle()
    gROOT.SetBatch(True)
    gStyle.SetPalette(55)

    print 'Overlaying %d average MinBias events'%opt.pu

    #later can add other options
    runOccupancyAnalysis(url=opt.input,avgPU=opt.pu,sdType=opt.sd,cellSize=opt.cell)

if __name__ == "__main__":
    main()


