#!/usr/bin/env python

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
def runOccupancyAnalysis(accMap,url='HGCSimHitsAnalysis.root',avgPU=0,sdType=0,noiseScale=0.0,mipEn=60,treeName='hgcSimHitsAnalyzer/HGC') :

    noiseSigma={0:0.075*noiseScale,
                1:0.3  *noiseScale,
                2:0.675*noiseScale,
                3:1.2  *noiseScale}
    if sdType==1 :
        noiseSigma[0]=noiseSigma[1]
        noiseSigma[1]=noiseSigma[2]
        noiseSigma[2]=noiseSigma[3]
        noiseSigma[3]=16*noiseSigma[0]
        
    thresholds=[0.4,1.0,2.0,5.0,10.0,25.0]
    etaVsEnergyPerLayer={}
    etaVsEnergy={}
    for cell in [0,1,2,3]:
        etaVsEnergy[cell]=ROOT.TH2F('espec_cell%d'%cell,';Pseudo-rapidity;Energy / MIP',15,1.5,3.0,500,0,20)
        etaVsEnergy[cell].SetDirectory(0)

    #prepare output
    foutUrl='occprofiles_pu%d_sd%d_noise%3.1f.root'%(avgPU,sdType,noiseScale)
    fout=TFile.Open(foutUrl,'RECREATE')
    tupleVars="layer:eta:cell:n"
    for thr in thresholds : tupleVars=tupleVars+':n_%d'%(thr*10)
    outTuple = ROOT.TNtuple("occ","occ",tupleVars)

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
                sector=Events.hit_sec[idep]-1

                #add
                edep=Events.hit_edep[idep]*1e6/mipEn
                bin=Events.hit_bin[idep]

                if sector<0: continue
                accMap[layer].fill(edep,bin,sector)

                
        #integrate up to 4x4 cells to determine energy deposits
        #  _______________________________
        # | (0,3)   (1,3)   (2,3)   (3,3) |
        # |_______________________        |
        # | (0,2)   (1,2)   (2,2) | (3,2) |
        # |_______________        |       |
        # | (0,1)   (1,1) | (2,1) | (3,1) |
        # |_______        |       |       |
        # | (0,0) | (1,0) | (2,0) | (3,0) |
        # |_______________________________|
        #
        for layer in accMap:
            if layer<0: continue

            accMap[layer].accumulator[0].Draw('colz')
            raw_input()
            
            xbin=accMap[layer].accumulator[0].GetXaxis().FindBin(0)
            for ybin in xrange(1,accMap[layer].accumulator[0].GetYaxis().GetNbins()-4,4):

                for sector in accMap[layer].accumulator:
                    enAcc=[]

                    #accMap[layer].accumulator[sector].Draw('colz')
                    totalEp_pos,totalEm_pos,totalEp_neg,totalEm_neg=0,0,0,0
                    totalEta,totalArea=0,0
                    for xStep,yStep in [[0,0],
                                        [0,1],[1,0],[1,1],
                                        [0,2],[1,2],[2,0],[2,1],[2,2],
                                        [0,3],[1,3],[2,3],[3,0],[3,1],[3,2],[3,3]]:

                        #positive x
                        totalEp_pos=totalEp_pos+accMap[layer].accumulator[sector].GetBinContent(xbin+xStep,ybin+yStep)
                        totalEm_pos=totalEm_pos+accMap[-layer].accumulator[sector].GetBinContent(xbin+xStep,ybin+yStep)

                        #negative x
                        totalEp_neg=totalEp_neg+accMap[layer].accumulator[sector].GetBinContent(xbin-1-xStep,ybin+yStep)
                        totalEm_neg=totalEm_neg+accMap[-layer].accumulator[sector].GetBinContent(xbin-1-xStep,ybin+yStep)

                        #eta for the integrated cell
                        totalEta=totalEta+accMap[layer].etaMap.GetBinContent(xbin+xStep,ybin+yStep)

                        #total area integrated
                        totalArea=totalArea+accMap[layer].accumulator[sector].GetXaxis().GetBinWidth(xbin+xStep)

                        #save at the edges
                        if xStep==yStep: enAcc.append( [xStep,
                                                        totalEta/((xStep+1)*(yStep+1)),
                                                        totalArea,
                                                        [totalEm_pos,totalEp_pos,totalEm_neg,totalEp_neg]] )

                    #fill the energy spectra histograms
                    print ybin,sector,enAcc[0][3]
                    for accStep in enAcc:
                        cell=accStep[0]
                        eta=accStep[1]
                        for en in accStep[3]:
                            noise=ROOT.gRandom.Gaus(0,noiseSigma[cell])
                            if noise<0: noise=0
                            etaVsEnergy[cell].Fill(eta,en+noise)

            for cell in etaVsEnergy:
                #etaVsEnergy[cell].Draw('colz')
                if not (layer in etaVsEnergyPerLayer):
                    etaVsEnergyPerLayer[layer]={}
                if not (cell in etaVsEnergyPerLayer[layer]):
                    etaVsEnergyPerLayer[layer][cell]=etaVsEnergy[cell].Clone('espec_layer%d_cell%d'%(layer,cell))
                    etaVsEnergyPerLayer[layer][cell].SetDirectory(0)
                etaVsEnergyPerLayer[layer][cell].Add(etaVsEnergy[cell])

                #now project for each eta and save in a smaller ntuple
                for etaBin in xrange(1,etaVsEnergy[cell].GetXaxis().GetNbins()+1):
                    pyH=etaVsEnergy[cell].ProjectionY("py",etaBin,etaBin)
                    n=pyH.Integral()
                    eta=etaVsEnergy[cell].GetXaxis().GetBinCenter(etaBin)
                    if eta>3 : continue
                    varsToFill=[layer,eta,cell,n]
                    for thr in thresholds:
                        startBin=pyH.GetXaxis().FindBin(thr)
                        varsToFill.append( pyH.Integral(startBin,pyH.GetXaxis().GetNbins()+1) )
                    outTuple.Fill(array("f",varsToFill))

            #all done with this layer
            accMap[layer].reset()
                 
    #all done with the inputs
    fin.Close()
    
    #save profiles obtained as function of y and y->eta map
    fout.cd()
    for layer in etaVsEnergyPerLayer:
        for cell in etaVsEnergyPerLayer[layer]:
            etaVsEnergyPerLayer[layer][cell].SetDirectory(fout)
            etaVsEnergyPerLayer[layer][cell].Write()
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


