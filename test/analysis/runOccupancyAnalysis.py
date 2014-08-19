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
loops over the events and collects the electron energies
"""
def runOccupancyAnalysis(accMap,fInUrl,tag,avgPU=0,sdType=0,mipEn=54.8,treeName='hgcSimHitsAnalyzer/HGC') :

    dataVolumePerLayer={}

    #prepare output
    foutUrl='%s_occ_pu%d_sd%d.root'%(tag,avgPU,sdType)
    fout=TFile.Open(foutUrl,'RECREATE')
    outTuple = ROOT.TNtuple("occ","occ","event:layer:cell:area:eta:e_sr0:e_sr1:e_sr2")
                            
    #analyze generated minimum bias events
    Events=ROOT.TChain(treeName)
    for f in fInUrl:
        Events.Add(f)
    print 'Found %d events in %d files'%(Events.GetEntries(),len(fInUrl))
    
    #generate pileup sets
    evRanges=generateInTimePUSets(Events.GetEntries(),avgPU)

    #useful to integrate over different cells
    cellInteg=CellIntegrator()
    
    #iterate over events
    iSetCtr=0
    for iSet in evRanges:
        if iSetCtr>1: break
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

                if layer%3==0: continue
                layer=layer/3

                #add
                edep=Events.hit_edep[idep]*1e6 #/mipEn
                bin=Events.hit_bin[idep]
                
                if sector<0: continue
                try:
                    accMap[layer].fill(edep,bin,sector)
                except:
                    pass

        #check what has been accumulated
        for layer in accMap:

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
                for cell in [1,2,3]:

                    centerX    = 0
                    centerXbin = accMap[layer].accumulator[sector].GetXaxis().FindBin(centerX)+cellInteg.center[cell][0]-1
                    centerX    = accMap[layer].accumulator[sector].GetXaxis().GetBinCenter(centerXbin)

                    startY     = accMap[layer].accumulator[sector].GetYaxis().GetXmin()
                    startYbin  = accMap[layer].accumulator[sector].GetYaxis().FindBin(startY)+cellInteg.center[cell][1]-1
                    startY     = accMap[layer].accumulator[sector].GetYaxis().GetBinCenter(startYbin)

                    for centerYbin in xrange(startYbin,accMap[layer].accumulator[sector].GetYaxis().GetNbins()-cellInteg.step[cell],cellInteg.step[cell]):
                        
                        nCellsInteg = {0:0, 1:0, 2:0 }
                        totalInCell = {0:0.,1:0.,2:0.}
                        etaInCell   = {0:0.,1:0.,2:0.}
                        for iXbin in xrange(cellInteg.integRange[cell][0],cellInteg.integRange[cell][1]+1):
                            for iYbin in xrange(cellInteg.integRange[cell][0],cellInteg.integRange[cell][1]+1):
                                srNumber=cellInteg.getSRNumber(iXbin,iYbin,cell)
                                nCellsInteg[srNumber] = nCellsInteg[srNumber] + 1
                                totalInCell[srNumber] = totalInCell[srNumber] + accMap[layer].accumulator[sector].GetBinContent(centerXbin+iXbin,centerYbin+iYbin)
                                etaInCell[srNumber]   = etaInCell[srNumber]   + accMap[layer].etaMap[sector].GetBinContent(centerXbin+iXbin,centerYbin+iYbin)

                        #fill ntuple
                        varsToFill=[iSetCtr,layer,cell,TMath.Power(cell*binWidth,2),etaInCell[0]/nCellsInteg[0],totalInCell[0],totalInCell[1],totalInCell[2]]
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
    parser.add_option('-i', '--in'         ,    dest='input'              , help='Input directory',                     default=None           )
    parser.add_option('-t', '--tag'        ,    dest='tag'                , help='Simulation tag',                      default="", )
    parser.add_option('-p', '--pu'         ,    dest='pu'                 , help='Average PU to overlay',               default=0,    type=int )
    parser.add_option('-s', '--sd'         ,    dest='sd'                 , help='Sensitive detector to analyse',       default=0,    type=int )    
    (opt, args) = parser.parse_args()

    if opt.input is None :
        parser.print_help()
        sys.exit(1)

    #plot formatting
    customROOTstyle()
    #gROOT.SetBatch(True)
    gROOT.SetBatch(False)
    gStyle.SetPalette(55)

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
    
    accMap=readSectorHistogramsFrom(fInUrl=fInUrl[0],sd=opt.sd)
    runOccupancyAnalysis(accMap,fInUrl=fInUrl,tag=opt.tag,avgPU=opt.pu,sdType=opt.sd)

if __name__ == "__main__":
    main()


