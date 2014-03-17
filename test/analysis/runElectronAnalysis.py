#!/usr/bin/env python

import ROOT
import sys
import optparse
import commands
import numpy as np
from array import array
from ROOT import *
from UserCode.HGCanalysis.PlotUtils import *
from UserCode.HGCanalysis.ElectronCandidate import *
from PileupUtils import *


"""
estimate transversed thickness correction for a given eta
"""
def getThicknessCorrectionForEta(etaVal):
    return 1./TMath.Cos(2*TMath.ATan(TMath.Exp(-etaVal)))

    
"""
loops over the events and collects the electron energies
"""
def runElectronAnalysis(url='particlegun.root',mbUrl='minbias.root',nPU=0,mipEn=54.8,treeName='hgcSimHitsAnalyzer/HGC') :

    customROOTstyle()
    #gROOT.SetBatch(True)
    gROOT.SetBatch(False)
    gStyle.SetPalette(56)

    #prepare the output
    fout=TFile.Open('ElectronAnalysis.root','RECREATE')
    fout.cd()
    output_tuple = TNtuple("etuple","etuple","en:pt:eta:sumEn:sumWEn:fitEn")
    output_tuple.SetDirectory(fout)
    
    #eta-phi display for energy deposits
    etaArray=array('d', sorted([-5]+np.arange(-3.4, -1.2, 0.05,dtype=np.double).tolist()+[-1.5,0,1.5]+np.arange(1.2,3.4, 0.05,dtype=np.double).tolist()+[5]) )
    etaPhiDisplay=TH2F('etaphi',';Pseudo-rapidity;#phi [rad];Energy',len(etaArray)-1,etaArray,128,-3.2,3.2)

    #open the minimum bias to overlay
    mbfin=TFile.Open(mbUrl)
    try:
        mbEvents=mbfin.Get(treeName)
        mbEvRanges=generateInTimePUSets(mbEvents.GetEntriesFast(),nPU)
    except:
        nPU=0
        print 'No minbias events will be mixed'

    #analyze generated events
    fin=TFile.Open(url)
    Events=fin.Get(treeName)
    eeBottomHalfWidth = fin.Get('hgcSimHitsAnalyzer/eeBottomHalfWidth')
    #baseOverburden=1*[0.21]+10*[0.5]+10*[0.8]+10*[1.2]
    for iev in xrange(6,Events.GetEntriesFast()) :

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
        genEta=Events.gen_eta[0]
        genPhi=Events.gen_phi[0]
        #thicknessSF=getThicknessCorrectionForEta(genEta)
        #localOverburden=[x0*thicknessSF for x0 in baseOverburden]
        localOverburden=100*[1.0]
        
        #get energy deposits info for ECAL
        if Events.nee==0 : continue
        edeps=len(localOverburden)*[0]
        edeps_etaphi={}
        for idep in xrange(0,Events.nee):
            edep=Events.ee_edep[idep]*1e6/mipEn
            #if edep<1 : continue
            layer=Events.ee_layer[idep]
            edeps[ layer-1 ] += edep
        
            #check in eta-phi
            gx=Events.ee_gx[idep]
            gy=Events.ee_gy[idep]
            gz=Events.ee_gz[idep]
            radius=TMath.Sqrt(gx*gx+gy*gy+gz*gz)
            p4=TLorentzVector(gx,gy,gz,radius)
            phi=p4.Phi()
            eta=p4.Eta()
            if idep==0 :
                print phi,Events.ee_x[idep],Events.ee_y[idep],Events.ee_subsec[idep],gy,gx,ROOT.TMath.ATan2(gy,gx)
            etaPhiDisplay.Fill(eta,phi,edep)

            #select around the electron a eta-phi square
            deta=eta-genEta
            dphi=TVector2.Phi_mpi_pi(phi-genPhi)
            #if TMath.Abs(deta)>1.0 or TMath.Abs(dphi)>1.5 : continue
            if not layer in edeps_etaphi: edeps_etaphi[layer]=[]
            edeps_etaphi[layer].append([deta,dphi,edep])


        #overlay the pileup
        if nPU>0:
            idxToOverlay=0 #int(ROOT.gRandom.Uniform(0,len(mbEvRanges)))
            
            print 'overlaying ',mbEvRanges[idxToOverlay][0],mbEvRanges[idxToOverlay][1]+1
            for imbev in xrange(mbEvRanges[idxToOverlay][0],mbEvRanges[idxToOverlay][1]+1):
                mbEvents.GetEntry(imbev)
                for idep in xrange(0,mbEvents.nee):
                    edep=mbEvents.ee_edep[idep]*1e6/mipEn
                    layer=mbEvents.ee_layer[idep]
                    edeps[ layer-1 ] += edep
                    
                    gx=mbEvents.ee_gx[idep]
                    gy=mbEvents.ee_gy[idep]
                    gz=mbEvents.ee_gz[idep]
                    radius=TMath.Sqrt(gx*gx+gy*gy+gz*gz)
                    p4=TLorentzVector(gx,gy,gz,radius)
                    phi=p4.Phi()
                    eta=p4.Eta()
                    etaPhiDisplay.Fill(eta,phi,edep) #*localOverburden[mbEvents.ee_layer[idep]-1]/localOverburden[0])
                    
                    #select around the electron a eta-phi square
                    deta=eta-genEta
                    dphi=TVector2.Phi_mpi_pi(phi-genPhi)
                    #if TMath.Abs(deta)>1.0 or TMath.Abs(dphi)>1.5 : continue
                    if not layer in edeps_etaphi:
                        edeps_etaphi[layer]=[]
                    edeps_etaphi[layer].append([deta,dphi,edep])


        #print edeps_etaphi
        #build the electron candidate
        ele=ElectronCandidate()
        ele.setGeneratorLevelInfo(genEn,genPt,genEta)
        ele.setLongitudinalProfile(edeps,localOverburden)
        ele.buildLongitudinalProfile(iev>10 and iev<20)
        ele.setTransverseProfile(edeps_etaphi)
        ele.buildTransverseProfile(iev<10)
        
        if iev<10:
            c=ROOT.TCanvas('c','c',500,500)
            c.SetTopMargin(0.05)
            etaPhiDisplay.Draw('lego2fb 0')

            phietaLine=ROOT.TPolyLine3D(3)
            phietaLine.SetPoint(0,genEta,-3.2,0)
            phietaLine.SetPoint(1,genEta,genPhi,0)
            phietaLine.SetPoint(2,-5,genPhi,0)
            phietaLine.SetLineStyle(9)
            phietaLine.SetLineWidth(2)
            phietaLine.Draw('same')

            #c.SetRightMargin(0.15)
            #etaPhiDisplay.Draw('colz')
            #etaPhiDisplay.GetZaxis().SetRangeUser(1,2e4)
            etaPhiDisplay.GetZaxis().SetLabelSize(0.03)
            etaPhiDisplay.GetZaxis().SetTitleSize(0.04)
            etaPhiDisplay.GetYaxis().SetLabelSize(0.03)
            etaPhiDisplay.GetYaxis().SetTitleSize(0.04)
            etaPhiDisplay.GetXaxis().SetLabelSize(0.03)
            etaPhiDisplay.GetXaxis().SetTitleSize(0.04)
            etaPhiDisplay.GetXaxis().SetTitleOffset(1.5)
            etaPhiDisplay.GetYaxis().SetTitleOffset(1.5)
            etaPhiDisplay.GetZaxis().SetTitleOffset(1.8)
            MyPaveText('CMS simulation, <PU>=%d'%nPU)
            pt=MyPaveText('E^{gen}=%3.0f GeV,p_{T}^{gen}=%3.0f GeV, #eta^{gen}=%3.1f #phi^{gen}=%3.1f'%(genEn,genPt,genEta,genPhi),0.1,0.9,0.8,0.92)
            pt.SetTextFont(42)
            pt.SetTextSize(0.03)
            c.Modified()
            c.Update()
            c.SaveAs('cmssw_%d_pu%d.png'%(iev,nPU))
            raw_input()
        etaPhiDisplay.Reset('ICE')
       

        #save output to tree
        values=[genEn,genPt,genEta,ele.totalRawEn,ele.totalEn,ele.totalFitEn]
        output_tuple.Fill(array("f",values))
        
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
    parser.add_option('-i', '--in'         ,    dest='input'              , help='Input file',              default=None)
    parser.add_option('-m', '--mb'         ,    dest='mbInput'            , help='Minimum bias input file', default='minbias.root')
    parser.add_option('-p', '--pu'         ,    dest='nPU'                , help='Average in-time PU to overlay', default=0, type=int)
    (opt, args) = parser.parse_args()

    if opt.input is None :
        parser.print_help()
        sys.exit(1)

    #later can add other options
    runElectronAnalysis(opt.input,opt.mbInput,opt.nPU)

if __name__ == "__main__":
    main()


