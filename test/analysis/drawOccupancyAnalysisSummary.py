#!/usr/bin/env python

import os
import sys
import optparse
import commands
from array import array
import numpy as np

from UserCode.HGCanalysis.PlotUtils import *

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
Return the weighted average and standard deviation.values, weights -- Numpy ndarrays with the same shape.
"""
def weighted_avg_and_std(values, weights):
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return average, TMath.Sqrt(variance)

"""
checks the input arguments and steers the analysis
"""
def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in'       ,    dest='input'          , help='Input file',  default=None)
    parser.add_option('-o', '--out'      ,    dest='output'         , help='Output dir', default=None)
    parser.add_option('-m', '--mip'      ,    dest='mip'            , help='MIP energy [keV]', default=55.1)
    (opt, args) = parser.parse_args()

    if opt.input is None or opt.output is None:
        parser.print_help()
        sys.exit(1)

    #plot formatting
    customROOTstyle()
    #gROOT.SetBatch(True)
    gROOT.SetBatch(False)
    gStyle.SetPalette(53)

    #prepare the output
    os.system('mkdir -p %s'%opt.output)

    fIn=TFile.Open(opt.input)
    occT=fIn.Get('occ')

    #get base cell area (mm^2)
    occT.GetEntry(0)
    baseArea=occT.area

    #get how many layers the sensitive detector has
    layersH=TH1F('layersH',';Layer #;',100,0,100)
    occT.Draw('TMath::Abs(layer) >> layersH','','goff')
    layers=[]
    for xbin in xrange(1,layersH.GetXaxis().GetNbins()+1):
        if layersH.GetBinContent(xbin)==0 : continue
        layers.append(xbin-1)
    layersH.Delete()


    cells=[1,2,3]
    thresholds=[0.4,1.0,5.0,10.,25.0]
    etaBins=[[1,2],[3,4],[5,6],[7,8],[9,10],[11,12],[13,15]]

    #prepare plots
    eprofH=TH2F('eprofH',';Pseudo-rapidity;(E/MIP) / coth(#eta);Hits',15,1.5,3.0,130,0,26)
    eprofH.GetZaxis().SetTitleOffset(-0.3)
    eprofH.GetZaxis().SetLabelSize(0.035)
    eprofH.Sumw2()
    eprofGr=TGraphErrors()
    eprofGr.SetMarkerStyle(20)
    eprofGr.SetMarkerColor(38)
    eprofGr.SetLineColor(38)
    eprofGr.SetLineWidth(2)
    occProfiles={}
    dataVolProfiles={}
    for thr in thresholds:
        nThrGr=len(occProfiles)
        occProfiles[thr]={}
        for cell in cells:
            occProfiles[thr][cell]={}
            if not( cell in dataVolProfiles) : dataVolProfiles[cell]={}
            for iEtaBin in xrange(0,len(etaBins)):
                nEtaGr=len(occProfiles[thr][cell])
                occProfiles[thr][cell][iEtaBin]=TGraphErrors()
                occProfiles[thr][cell][iEtaBin].SetMarkerStyle(20+nEtaGr)
                occProfiles[thr][cell][iEtaBin].SetMarkerColor(1+nEtaGr%3)
                occProfiles[thr][cell][iEtaBin].SetLineColor(1+nEtaGr%3)
                occProfiles[thr][cell][iEtaBin].SetLineWidth(2)
                occProfiles[thr][cell][iEtaBin].SetTitle('%3.1f-%3.1f'%(eprofH.GetXaxis().GetBinLowEdge(etaBins[iEtaBin][0]),
                                                                        eprofH.GetXaxis().GetBinUpEdge(etaBins[iEtaBin][1])))
                if not(iEtaBin in dataVolProfiles[cell]) : dataVolProfiles[cell][iEtaBin]=occProfiles[thr][cell][iEtaBin].Clone()
                
    #base canvas
    c=TCanvas('c','c',500,500)

    #iterate over layers
    for layer in layers:
        for cell in cells:

            #energy profile versus eta
            eprofH.Reset('ICE')
            occT.Draw('(e_sr0/%f)*TMath::TanH(eta):eta >> eprofH'%(opt.mip),'layer==%d && cell==%d'%(layer,cell),'goff')
            eprofH.Draw('colz')
            MyPaveText('CMS simulation, <PU>=140 (14 TeV)')
            ptxt=MyPaveText('#bf{[Layer #%d]} cell area=%3.1f mm^{2}'%(layer,baseArea*cell),0.12,0.9,0.9,0.93)
            ptxt.SetTextSize(0.035)
            ptxt.SetFillColor(0)
            ptxt.SetFillStyle(3001)
            
            #project in eta slices
            eprofGr.Set(0)
            eDistsH=[]
            for iEtaBin in xrange(0,len(etaBins)):

                #coordinates
                etaMin=eprofH.GetXaxis().GetBinLowEdge(etaBins[iEtaBin][0])
                etaMax=eprofH.GetXaxis().GetBinUpEdge(etaBins[iEtaBin][1])
                ieta=0.5*(etaMax+etaMin)
                deta=etaMax-etaMin

                #project
                eH=eprofH.ProjectionY('eH',etaBins[iEtaBin][0],etaBins[iEtaBin][1])
                eH.SetTitle('%3.1f-%3.1f'%(etaMin,etaMax))
                fixExtremities(eH)
                totalHits=eH.Integral()
                eDistsH.append( eH.Clone('eprof_%d'%iEtaBin) )
                if totalHits==0 : continue

                #get average deposit
                meanEn,meanEnErr=eH.GetMean(),eH.GetMeanError()
                np=eprofGr.GetN()
                eprofGr.SetPoint(np,ieta,meanEn)
                eprofGr.SetPointError(np,0.5*deta,meanEnErr)
                
                #compute occupancies
                for thr in thresholds:
                    ybin=eH.GetXaxis().FindBin(thr)
                    cellOccErr=ROOT.Double(0)
                    cellOcc=eH.IntegralAndError(ybin,eH.GetXaxis().GetNbins(),cellOccErr)
                    np=occProfiles[thr][cell][iEtaBin].GetN()
                    occProfiles[thr][cell][iEtaBin].SetPoint(np,layer,cellOcc*100/totalHits)
                    occProfiles[thr][cell][iEtaBin].SetPointError(np,0,cellOccErr*100/totalHits)
                    
                #data volume
                dataSent=[]
                evCount=[]
                for ybin in xrange(1,eH.GetXaxis().GetNbins()+1):
                    evCount.append( eH.GetBinContent(ybin) )
                    dataSent.append( dataBitsToSend(ieta, eH.GetXaxis().GetBinCenter(ybin))[0] )
                np=dataVolProfiles[cell][iEtaBin].GetN()
                avgDataVol,stdDataVol=weighted_avg_and_std(dataSent,evCount)
                dataVolProfiles[cell][iEtaBin].SetPoint(np,layer,avgDataVol)
                dataVolProfiles[cell][iEtaBin].SetPointError(np,0,stdDataVol/TMath.Sqrt(totalHits))

                #all done here
                eH.Delete()
                
            #save energy profile canvas
            eprofGr.Draw('p')
            c.SetRightMargin(0.12)
            c.SetLogz(True)
            c.SetLogy(False)
            c.Modified()
            c.Update()
            c.SaveAs('%s/eprof2D_layer%d_cell%d.png'%(opt.output,layer,cell))

            #energy distributions at a given #eta
            c.Clear()
            c.SetLogz(False)
            c.SetRightMargin(0.05)
            
            drawOpt='hist'
            leg=TLegend(0.15,0.7,0.45,0.88)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.SetTextFont(42)
            leg.SetTextSize(0.035)
            leg.SetHeader('Pseudo-rapidity')
            for ih in xrange(len(eDistsH)-1,0,-2):
                eDistsH[ih].Draw(drawOpt)
                drawOpt='hist same'
                eDistsH[ih].GetYaxis().SetTitle('Hits')
                ctr=(len(eDistsH)-ih)/2
                eDistsH[ih].SetLineColor(1+ctr/2)
                eDistsH[ih].SetLineStyle(1+ctr%2)
                eDistsH[ih].SetLineWidth(2)
                leg.AddEntry(eDistsH[ih],eDistsH[ih].GetTitle(),'l')

            leg.Draw()
            MyPaveText('CMS simulation, <PU>=140 (14 TeV)')
            ptxt=MyPaveText('#bf{[Layer #%d]} cell area=%3.1f mm^{2}'%(layer,baseArea*cell),0.12,0.9,0.9,0.93)
            ptxt.SetTextSize(0.035)
            ptxt.SetFillColor(0)
            ptxt.SetFillStyle(3001)
            
            #zoom in
            c.cd()
            subc=TPad('zoom','zoom',0.47,0.45,0.92,0.9)
            subc.Draw()
            subc.cd();
            subc.SetFillStyle(0)
            subc.SetBottomMargin(0.25)
            subc.SetLeftMargin(0.25)
            drawOpt='hist'
            zoomDistsH=[]
            for ih in xrange(len(eDistsH)-1,0,-2):
                nzoom=len(zoomDistsH)
                zoomDistsH.append( eDistsH[ih].Clone('%s_zoom'%eDistsH[ih].GetName()) )
                zoomDistsH[ nzoom ].Draw(drawOpt)
                zoomDistsH[ nzoom ].GetXaxis().SetRangeUser(0.4,3)
                drawOpt='hist same'

            c.cd()
            c.SetLogy()
            c.Modified()
            c.Update()
            c.SaveAs('%s/eprof1D_layer%d_cell%d.png'%(opt.output,layer,cell))

    #occupancy summaries
    for thr in occProfiles:
        for cell in occProfiles[thr]:

            c.Clear()
            c.SetRightMargin(0.05)
            c.SetLogz(False)
            c.SetLogy(False)
            drawOpt='acp'
            leg=TLegend(0.8,0.72,0.95,0.93)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.SetTextFont(42)
            leg.SetTextSize(0.035)
            leg.SetHeader('Pseudo-rapidity')
            for eta in occProfiles[thr][cell]:
                #if occProfiles[thr][cell][eta].GetN()<layers : continue
                occProfiles[thr][cell][eta].Draw(drawOpt)
                occProfiles[thr][cell][eta].GetXaxis().SetTitle("Layer")
                occProfiles[thr][cell][eta].GetYaxis().SetTitle("Occupancy (%)")
                occProfiles[thr][cell][eta].GetYaxis().SetRangeUser(0,100)
                drawOpt='cp'
                leg.AddEntry(occProfiles[thr][cell][eta],occProfiles[thr][cell][eta].GetTitle(),'cp')

            leg.Draw()
            MyPaveText('CMS simulation, <PU>=140 (14 TeV)')
            MyPaveText('Cell area=%3.1f mm^{2}, Threshold=%3.1f/MIP'%(baseArea*cell,thr),0.12,0.9,0.9,0.93).SetTextSize(0.035)

            c.Modified()
            c.Update()
            c.SaveAs('%s/occ_thr%d_cell%d.png'%(opt.output,10*thr,cell))
            c.SetLogy(True)
            for eta in occProfiles[thr][cell]: occProfiles[thr][cell][eta].GetYaxis().SetRangeUser(0,100)
            c.Modified()
            c.Update()
            c.SaveAs('%s/occ_thr%d_cell%d_log.png'%(opt.output,10*thr,cell))

    #data volume summaries
    for cell in dataVolProfiles:
        c.Clear()
        c.SetRightMargin(0.05)
        c.SetLogz(False)
        c.SetLogy(False)
        drawOpt='acp'
        leg=TLegend(0.8,0.72,0.95,0.93)
        leg.SetHeader('Pseudo-rapidity')
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        for eta in dataVolProfiles[cell]:
            dataVolProfiles[cell][eta].Draw(drawOpt)
            dataVolProfiles[cell][eta].GetXaxis().SetTitle("Layer")
            dataVolProfiles[cell][eta].GetYaxis().SetTitle("Average data sent / cell (bit)")
            dataVolProfiles[cell][eta].GetYaxis().SetRangeUser(0,3)
            drawOpt='cp'
            leg.AddEntry(dataVolProfiles[cell][eta],dataVolProfiles[cell][eta].GetTitle(),'cp')

        leg.Draw()
        MyPaveText('CMS simulation, <PU>=140 (14 TeV)')
        MyPaveText('Cell area=%3.1f mm^{2}'%(baseArea*cell),0.12,0.9,0.9,0.93).SetTextSize(0.035)

        c.Modified()
        c.Update()
        c.SaveAs('%s/datavol_cell%d.png'%(opt.output,cell))

            
    fIn.Close()
    

if __name__ == "__main__":
    main()


