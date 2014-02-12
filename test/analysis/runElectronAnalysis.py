#!/usr/bin/env python

import ROOT
import sys
import optparse
import commands
import array
from ROOT import *
from UserCode.HGCanalysis.PlotUtils import *
from UserCode.HGCanalysis.ElectronCandidate import *



"""
estimate transversed thickness correction for a given eta
"""
def getThicknessCorrectionForEta(etaVal):
    return 1./TMath.Cos(2*TMath.ATan(TMath.Exp(-etaVal)))

    
"""
loops over the events and collects the electron energies
"""
def runElectronAnalysis(url='HGCSimHitsAnalysis.root',mipEn=54.6,treeName='hgcSimHitsAnalyzer/HGC') :

    customROOTstyle()
    gROOT.SetBatch(True)
    gStyle.SetPalette(55)

    #prepare the output
    fout=TFile.Open('ElectronAnalyis.root','RECREATE')
    fout.cd()
    output_tuple = TNtuple("etuple","etuple","en:pt:eta:sumEn:sumWEn:fitEn")

    #analyze generated events
    fin=TFile.Open(url)
    Events=fin.Get(treeName)
    baseOverburden=10*[0.5]+10*[0.8]+10*[1.2]
    for iev in xrange(0,Events.GetEntriesFast()) :

        #some printout to know where we are
        if iev%100 == 0 :
            sys.stdout.write("\r[%d/100]"%(100.*iev/Events.GetEntriesFast()))
            sys.stdout.flush()
            
        Events.GetEntry(iev)

        #get generator level info
        if Events.ngen==0: continue
        id=Events.gen_id[0]
        if abs(id) != 11 : continue
        genEn=Events.gen_en[0]
        genPt=Events.gen_pt[0]
        genEta=abs(Events.gen_eta[0])
        thicknessSF=getThicknessCorrectionForEta(genEta)
        localOverburden=[x0*thicknessSF for x0 in baseOverburden]

        #get energy deposits info for ECAL
        if Events.nee==0 : continue
        edeps=len(localOverburden)*[0]
        for idep in xrange(0,Events.nee):
            edep=Events.ee_edep[idep]*1e6/mipEn
            if edep<1 : continue
            edeps[ Events.ee_layer[idep]-1 ] += edep

        #build the electron candidate
        ele=ElectronCandidate()
        ele.setGeneratorLevelInfo(genEn,genPt,genEta)
        ele.setLongitudinalProfile(edeps,localOverburden)
        ele.buildLongitudinalProfile(iev<10)

        #save output to tree
        values=[genEn,genPt,genEta,ele.totalRawEn,ele.totalEn,ele.totalFitEn]
        output_tuple.Fill(array.array("f",values))
        
    fout.cd()
    output_tuple.Write()
    fout.Close()        

    fin.Close()
    



"""
checks the input arguments and steers the analysis
"""
def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in'         ,    dest='input'              , help='Input file', default=None)
    (opt, args) = parser.parse_args()

    if opt.input is None :
        parser.print_help()
        sys.exit(1)

    #later can add other options
    runElectronAnalysis(opt.input)

if __name__ == "__main__":
    main()


