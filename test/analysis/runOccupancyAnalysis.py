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
    h.Sumw2()
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

            layer=(layer-1)/sdCtrPeriod+1

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
def runOccupancyAnalysis(url='HGCSimHitsAnalysis.root',avgPU=0,sdType=0,cellSize=10,iniSetCtr=0,finalSetCtr=-1,mipEn=54.8,treeName='hgcSimHitsAnalyzer/HGC') :

    #get templates 
    




    #analyze generated minimum bias events
    fin=TFile.Open(url)
    Events=fin.Get(treeName)
    


    

    rangePostFix=''
    if finalSetCtr>iniSetCtr : rangePostFix='_%dto%d'%(iniSetCtr,finalSetCtr)
    foutUrl='occprofiles_pu%d_sd%d_cell%d%s.root'%(avgPU,sdType,cellSize,rangePostFix)

    fout=TFile.Open(foutUrl,'RECREATE')
    outTuple = ROOT.TNtuple("occ","occ","layer:eta:en:iso1:iso2:iso3")
    
    #base templates in global coordinates and in local coordinates
    sdCtrPeriod=1
    subLayerToSample=[0]
    if sdType in [0,1]:
        sdCtrPeriod=3
        subLayerToSample=[0,1]
    layerHistos=getLayerSubSectorHistos(fin,cellSize,sdType,sdCtrPeriod)

    #
    # histograms for analysis
    #
    layerEnAverage={}
    layerOccAverage={}
    for l in layerHistos:

        #energy distributions
        ybins=layerHistos[l][0].GetYaxis().GetNbins()
        ymin=layerHistos[l][0].GetYaxis().GetXmin()
        ymax=layerHistos[l][0].GetYaxis().GetXmax()
        layerEnAverage[l]=TH2F('avgen_layer%d'%(l),';<Energy / MIP>; y [mm]; Events',100,0.,5.,ybins,ymin,ymax)
        layerEnAverage[l].SetDirectory(0)
        layerEnAverage[l].Sumw2()
            
        #occupancy profiles
        for thr in [0.4,1,5,10,25,50]:
            if not thr in layerOccAverage : layerOccAverage[thr]={}
            layerOccAverage[thr][l]=layerHistos[l][0].Clone('avgocc_%dmip_layer%d'%(thr*10,l))
            layerOccAverage[thr][l].SetDirectory(0)

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
        lbin=sdTranslX.GetXaxis().FindBin((l-1)*sdCtrPeriod+1)
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
    iGunCtr=0
    for iSet in evRanges:
        
        iSetCtr=iSetCtr+1
        sys.stdout.write( '\r Status [%d/%d]'%(iSetCtr,len(evRanges)) )
        sys.stdout.flush()

        if iSetCtr<iniSetCtr: continue
        if finalSetCtr>iniSetCtr :
            if iSetCtr>finalSetCtr: continue

                
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

                #print layer
                subLayer=(layer-1)%sdCtrPeriod
                layer=(layer-1)/sdCtrPeriod+1

                #fill base histos                
                if subLayer in subLayerToSample:
                    layerHistos[ layer ][ sector ].Fill(x,y,edep)
                
        #sum the histos over all the sectors
        for l in layerHistos:

            isec=0
            for secH in layerHistos[l]:
                isec=isec+1
                centerXbin=secH.GetXaxis().FindBin(0)                
                maxXbins=secH.GetXaxis().GetNbins()
                maxYbins=secH.GetYaxis().GetNbins()
            
                for ybin in xrange(1,maxYbins+1) :
                    
                    yval=secH.GetYaxis().GetBinCenter(ybin)

                    for xbin in xrange(1,maxXbins+1) :
                        
                        xval=secH.GetXaxis().GetBinCenter(xbin)
                        enInCell=secH.GetBinContent(xbin,ybin)
                        
                        #energy spectrum in single cell
                        if xbin == centerXbin:

                            #average energy, over sector
                            layerEnAverage[l].Fill(yval,enInCell)

                            #isolation
                            sumInRing={}
                            nCellsInSum={}
                            for ring in [1,2,3]:
                                sumInRing[ring]=0
                                nCellsInSum[ring]=0
                                for ix in xrange(xbin-ring,xbin+ring+1):
                                    for iy in xrange(ybin-ring,ybin+ring+1):
                                        if not (ix==xbin and iy==ybin) :
                                            if ix>0 and ix<=maxXbins and iy>0 and iy<=maxYbins : 
                                                sumInRing[ring]=sumInRing[ring]+secH.GetBinContent(ix,iy)
                                                nCellsInSum[ring]=nCellsInSum[ring]+1

                            etaBin=cellEtaMap[l].FindBin(yval)
                            eta=cellEtaMap[l].GetBinContent(etaBin)

                            valsToStore=[l,eta,enInCell,sumInRing[1],sumInRing[2]-sumInRing[1],sumInRing[3]-sumInRing[2]]
                            outTuple.Fill( array('f',valsToStore) )

                        #occupancy
                        for thr in layerOccAverage:
                            if enInCell<thr : continue
                            layerOccAverage[thr][l].Fill(xval,yval,1)
                 
    #all done with the inputs
    fin.Close()
    
    #save profiles obtained as function of y and y->eta map
    for l in layerHistos:
        cellEtaMap[l].SetDirectory(fout)
        cellEtaMap[l].Write()

        avgWeight=1.0/iSetCtr
        secWeight=1./len(layerHistos[l])
        finalWeight=avgWeight*secWeight

        layerEnAverage[l].Scale(finalWeight)
        layerEnAverage[l].SetDirectory(fout)
        layerEnAverage[l].Write()

        for thr in layerOccAverage:
            layerOccAverage[thr][l].Scale(finalWeight)
            layerOccAverage[thr][l].SetDirectory(fout)
            layerOccAverage[thr][l].Write()

    outTuple.Write()
    fout.Close()
    print 'Results stored in %s'%foutUrl


"""
checks the input arguments and steers the analysis
"""
def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in'         ,    dest='input'              , help='Input file',                          default=None           )
    parser.add_option('-p', '--pu'         ,    dest='pu'                 , help='Average PU to overlay',               default=0,    type=int )
    parser.add_option('-s', '--sd'         ,    dest='sd'                 , help='Sensitive detector to analyse',       default=0,    type=int )
    parser.add_option('-c', '--cell'       ,    dest='cell'               , help='cell size [mm]',                      default=10,   type=int )
    parser.add_option('-r', '--range'      ,    dest='range'              , help='first:last event to process',         default='0:-1')
    (opt, args) = parser.parse_args()

    if opt.input is None :
        parser.print_help()
        sys.exit(1)


    #plot formatting
    customROOTstyle()
    gROOT.SetBatch(True)
    #gROOT.SetBatch(False)
    gStyle.SetPalette(55)

    print 'Overlaying %d average MinBias events'%opt.pu
    iniSetCtr=int(opt.range.split(':')[0])
    finalSetCtr=int(opt.range.split(':')[1])
    print 'Analysing %d to %d events'%(iniSetCtr,finalSetCtr)
       
    #later can add other options
    runOccupancyAnalysis(url=opt.input,avgPU=opt.pu,sdType=opt.sd,cellSize=opt.cell,iniSetCtr=iniSetCtr,finalSetCtr=finalSetCtr)

if __name__ == "__main__":
    main()


