#!/usr/bin/env python

import os
import sys
import optparse
import commands
from array import array
import numpy as np

from UserCode.HGCanalysis.PlotUtils import *
from UserCode.HGCanalysis.PileupUtils import *
from UserCode.HGCanalysis.HGCLayerUtils import *

from ROOT import *

"""
Marcello's algorithm
"""
def dataBitsToSend(eta,nMips):
    if nMips==0 : return 0,0
    nBits=0
    region=0
    if eta<2.0:
        region=1
        if nMips<4:    nBits=1
        elif nMips<10: nBits=4+2
        else:          nBits=10
    elif eta<2.5:
        region=2
        if nMips<10:   nBits=1
        elif nMips<96: nBits=4+2
        else:          nBits=10
    else:
        region=3
        if nMips<25:    nBits=1
        elif nMips<192: nBits=4+2
        else:           nBits=10
    return nBits,region

"""
loops over the events and collects the electron energies
"""
def runOccupancyAnalysis(accMap,url='HGCSimHitsAnalysis.root',avgPU=0,sdType=0,noiseScale=0.0,mipEn=54.8,treeName='hgcSimHitsAnalyzer/HGC') :

    noiseSigma={0:0.075*noiseScale,
                1:0.3  *noiseScale,
                2:0.675*noiseScale,
                3:1.2  *noiseScale}
    if sdType==1 :
        noiseSigma[0]=noiseSigma[1]
        noiseSigma[1]=noiseSigma[2]
        noiseSigma[2]=noiseSigma[3]
        noiseSigma[3]=16*noiseSigma[0]
        
    
    dataVolumePerLayer={}

    #prepare output
    foutUrl=os.path.basename(url)
    foutUrl=foutUrl.replace('.root','_occ_pu%d_sd%d.root'%(avgPU,sdType))
    fout=TFile.Open(foutUrl,'RECREATE')
    outTuple = ROOT.TNtuple("occ","occ","event:layer:cell:eta:edep:edepRing")
                            
    #analyze generated minimum bias events
    fin=TFile.Open(url)
    Events=fin.Get(treeName)
    
    #generate pileup sets
    evRanges=generateInTimePUSets(Events.GetEntriesFast(),avgPU)
    
    #iterate over events
    iSetCtr=0
    for iSet in evRanges:
        
        iSetCtr=iSetCtr+1
        sys.stdout.write( '\r Status [%d/%d]'%(iSetCtr,len(evRanges)) )
        sys.stdout.flush()

        #sum up N events
        print iSet[0],iSet[1]
        for iev in xrange(iSet[0],iSet[1]+1):

            #get hits for new event
            Events.GetEntry(iev)
            if Events.nhits==0 : continue

            for idep in xrange(0,Events.nhits):

                #filter on sensitive detector
                if Events.hit_type[idep]!=sdType : continue

                #layer and sector info
                layer=Events.hit_layer[idep]
                sector=Events.hit_sec[idep]

                #add
                edep=Events.hit_edep[idep]*1e6 #/mipEn
                bin=Events.hit_bin[idep]

                if sector<0: continue
                accMap[layer].fill(edep,bin,sector)


        #check what has been accumulated
        for layer in accMap:

            #accMap[layer].accumulator[0].Draw('colz')
            #raw_input()

            #data volume estimate
            #for each cell compute the number of bits to send and integrate for a sector or specific eta's
            #if not (layer in dataVolumePerLayer):
            #    dataVolumePerLayer[layer]={}
            #    for reg in [0,1,2,3]:
            #       dataVolumePerLayer[layer][reg]=ROOT.TH1F('datavol_layer%d_reg_%d'%(layer,reg),';Data volume/bunch crossing (kb);Events x 20^{o} sector;',200,0,2)
            #        dataVolumePerLayer[layer][reg].SetDirectory(0)
            #for sector in accMap[layer].accumulator:
            #    totalBits={0:0,1:0,2:0,3:0}
            #    for ybin in xrange(1,accMap[layer].accumulator[sector].GetYaxis().GetNbins()+1):
            #        for xbin in xrange(1,accMap[layer].accumulator[sector].GetXaxis().GetNbins()+1):
            #            bitsToSend,ireg=dataBitsToSend( accMap[layer].etaMap[sector].GetBinContent(xbin,ybin), accMap[layer].accumulator[sector].GetBinContent(xbin,ybin) )
            #            if ireg==0 : continue
            #            totalBits[ireg]=totalBits[ireg]+bitsToSend
            #            totalBits[0]=totalBits[0]+bitsToSend
            #    for ireg in totalBits:
            #        dataVolumePerLayer[layer][ireg].Fill(float(totalBits[ireg])/1024.)

            #occupancy estimate
            #integrate up to 3x3 cells with increasing size (multiple of the base bin width)
            #and determine energy deposits and isolation
            #  _______________________ 
            # | (0,2)   (1,2)   (2,2) |
            # |       ________        |
            # | (0,1) | (1,1) | (2,1) |
            # |       ________|       |
            # | (0,0)   (1,0)   (2,0) | 
            # |_______________________|
            #
            #circulate over the sectors
            for sector in accMap[layer].accumulator :
                
                binWidth=accMap[layer].accumulator[sector].GetXaxis().GetBinWidth(1)

                #test different cell sizes
                for cellSize in [binWidth,2*binWidth,3*binWidth]:
                
                    halfStep=3*cellSize/2.
                    nCellsToInteg=int(halfStep/binWidth)+2
                
                    #probe at +/-x
                    for centerX in [-2*cellSize,2*cellSize]:
                        centerXbin=accMap[layer].accumulator[sector].GetXaxis().FindBin(centerX)
                        centerY=accMap[layer].accumulator[sector].GetYaxis().GetXmin()+2*cellSize

                        #go up in y
                        while centerY<accMap[layer].accumulator[sector].GetYaxis().GetXmax()-2*cellSize:

                            centerYbin=accMap[layer].accumulator[sector].GetYaxis().FindBin(centerY)

                            eta=0
                            nCellsInteg, nNeighbourCellsInteg=0,0
                            totalInCell,totalInNeighbours=0,0
                            for iXbin in xrange(-nCellsToInteg,nCellsToInteg):
                                for iYbin in xrange(-nCellsToInteg,nCellsToInteg):
                                    edep=accMap[layer].accumulator[sector].GetBinContent(centerXbin+iXbin,centerYbin+iYbin)
                                    ix=accMap[layer].accumulator[sector].GetXaxis().GetBinCenter(centerXbin+iXbin)
                                    iy=accMap[layer].accumulator[sector].GetYaxis().GetBinCenter(centerYbin+iYbin)

                                    #discard if outside the required range
                                    if TMath.Abs(ix-centerX)>halfStep or TMath.Abs(iy-centerY)>halfStep :
                                        continue
                                    elif TMath.Abs(ix-centerX)<cellSize/2 and TMath.Abs(iy-centerY)<cellSize/2 :
                                        totalInCell=totalInCell+edep
                                        nCellsInteg=nCellsInteg+1
                                        eta=eta+accMap[layer].etaMap[sector].GetBinContent(centerXbin+iXbin,centerYbin+iYbin)
                                    else:
                                        totalInNeighbours=totalInNeighbours+edep
                                        nNeighbourCellsInteg=nNeighbourCellsInteg+1
                                        
                            #move up
                            centerY=centerY+2*halfStep

                            #no noise for the moment
                            # noise=ROOT.gRandom.Gaus(0,noiseSigma[cell])

                            #if something has been integrated analyse it
                            if nCellsInteg==0: continue

                            #average eta of this cell
                            eta=eta/nCellsInteg

                            #fill ntuple
                            varsToFill=[iSetCtr,layer,cellSize/binWidth,eta,totalInCell,totalInNeighbours]
                            outTuple.Fill(array("f",varsToFill))                        

            #all done with this layer
            accMap[layer].reset()
                 
    #all done with the inputs
    fin.Close()
    
    #save profiles obtained as function of y and y->eta map
    fout.cd()
    #for layer in dataVolumePerLayer:
    #    for reg in dataVolumePerLayer[layer]:
    #        dataVolumePerLayer[layer][reg].SetDirectory(fout)
    #        dataVolumePerLayer[layer][reg].Write()        
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
    parser.add_option('-n', '--noise'      ,    dest='noiseScale'         , help='noise scale factor',                   default=0,   type=float )
    (opt, args) = parser.parse_args()

    if opt.input is None :
        parser.print_help()
        sys.exit(1)

    #plot formatting
    #customROOTstyle()
    #gROOT.SetBatch(True)
    #gROOT.SetBatch(False)
    #gStyle.SetPalette(55)
    
    accMap=readSectorHistogramsFrom(fInUrl=opt.input,sd=opt.sd)
    runOccupancyAnalysis(accMap,url=opt.input,avgPU=opt.pu,sdType=opt.sd,noiseScale=opt.noiseScale)

if __name__ == "__main__":
    main()


