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
    nBits=0
    if nMips==0 : return nBits
    if eta<2.0:
        if nMips<4:    nBits=1
        elif nMips<10: nBits=4+2
        else:          nBits=10
    elif eta<2.5:
        if nMips<10:   nBits=1
        elif nMips<96: nBits=4+2
        else:          nBits=10
    else:
        if nMips<25:    nBits=1
        elif nMips<192: nBits=4+2
        else:           nBits=10
    return nBits
    
"""
Return the weighted average and standard deviation.values, weights -- Numpy ndarrays with the same shape.
"""
def weighted_avg_and_std(values, weights):
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return average, TMath.Sqrt(variance)

"""
"""
def drawOccupancyAnalysisSummary(input,output,mip,addCut='',noiseLevel=0,pfix='') :

    #load utils
    gSystem.Load("libFWCoreFWLite.so");
    from ROOT import AutoLibraryLoader
    AutoLibraryLoader.enable()
    gSystem.Load("libUserCodeHGCanalysis.so")

    #prepare the output
    os.system('mkdir -p %s'%output)

    fIn=TFile.Open(input)
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

    #peak finder
    mipPeakFinderFunc=TF1('mipPeak','landau',0.4,1.6)
    
    cells=[1]
    #cells=[1,2,3]
    #thresholds=[0.4,1.0,5.0,10.,25.0]
    thresholds=[0.4]
    etaBins=[[1,3],[4,5],[6,7],[8,9],[10,11],[12,13],[14,15]]

    
    #generate code to compress data
    from ROOT import getHuffmanCodesFrom
    espectrumH=TH1F('espectrumH',';(E/MIP) / coth(#eta);Hits',15,0.4,6.4)
    espectrumH.Sumw2()
    huffmanCodes={}
    for cell in cells:
        espectrumH.Reset('ICE')
        occT.Draw('e_sr0/%f*TMath::TanH(eta) >> espectrumH'%mip,'cell==%d'%cell,'goff')
        huffmanCodes[cell]=getHuffmanCodesFrom(espectrumH)
        for bin in huffmanCodes[cell]:
            print bin,huffmanCodes[cell][bin].val,huffmanCodes[cell][bin].nbits
        
    #prepare plots
    eprofH=TH2F('eprofH',';Pseudo-rapidity;(E/MIP) / coth(#eta);Hits',15,1.5,2.8,130,0,26)
    eprofH.GetZaxis().SetTitleOffset(-0.3)
    eprofH.GetZaxis().SetLabelSize(0.035)
    eprofH.Sumw2()
    eprofGr=TGraphErrors()
    eprofGr.SetMarkerStyle(20)
    eprofGr.SetMarkerColor(38)
    eprofGr.SetLineColor(38)
    eprofGr.SetLineWidth(2)
    occProfiles={}
    dataVolH=TH2F('datavolH',';Pseudo-rapidity;#bits sent/bunch crossing;Hits',15,1.5,2.8,16,0,16)
    dataVolH.GetZaxis().SetTitleOffset(-0.3)
    dataVolH.GetZaxis().SetLabelSize(0.035)
    dataVolH.Sumw2()
    dataVolProfiles={}
    uncalibMIPprofiles={}
    for thr in thresholds:
        nThrGr=len(occProfiles)
        occProfiles[thr]={}
        for cell in cells:
            occProfiles[thr][cell]={}
            if not( cell in uncalibMIPprofiles):
                dataVolProfiles[cell]={}
                uncalibMIPprofiles[cell]={}
            for iEtaBin in xrange(0,len(etaBins)):
                nEtaGr=len(occProfiles[thr][cell])
                occProfiles[thr][cell][iEtaBin]=TGraphErrors()
                occProfiles[thr][cell][iEtaBin].SetMarkerStyle(20+nEtaGr)
                occProfiles[thr][cell][iEtaBin].SetMarkerColor(1+nEtaGr%3)
                occProfiles[thr][cell][iEtaBin].SetLineColor(1+nEtaGr%3)
                occProfiles[thr][cell][iEtaBin].SetLineWidth(2)
                occProfiles[thr][cell][iEtaBin].SetTitle('%3.1f-%3.1f'%(eprofH.GetXaxis().GetBinLowEdge(etaBins[iEtaBin][0]),
                                                                        eprofH.GetXaxis().GetBinUpEdge(etaBins[iEtaBin][1])))
                if not(iEtaBin in uncalibMIPprofiles[cell]):
                    dataVolProfiles[cell][iEtaBin]=occProfiles[thr][cell][iEtaBin].Clone()
                    uncalibMIPprofiles[cell][iEtaBin]=occProfiles[thr][cell][iEtaBin].Clone()
                
    #base canvas
    c=TCanvas('c','c',500,500)

    #iterate over layers
    for layer in layers:
        for cell in cells:

            #create an entry list with these cuts
            finalCut='layer==%d && cell==%d'%(layer,cell)
            if len(addCut) : finalCut=finalCut + ' && ' + addCut            
            occT.Draw('>>elist', finalCut,'entrylist')
            elist = ROOT.gDirectory.Get("elist")

            #fill the base histograms 
            datavolH.Reset('ICE')
            eprofH.Reset('ICE')
            for ientry in xrange(0,elist.GetN()):
                tentry=elist.GetEntry(ientry) 
                occT.GetEntry( tentry )
                corrEn=(occT.e_sr0/mip)*TMath.TanH(occT.eta)
                if noiseLevel>0 : corrEn = corrEn + gRandom.Gaus(0,noiseLevel)
                datavolH.Fill(TMath.Abs(occT.eta),dataBitsToSend(occT.eta,corrEn))
                eprofH.Fill(TMath.Abs(occT.eta),corrEn)

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
                
                #fit near the MIP peak
                mipPeakFinderFunc.SetParameter(1,1)
                eH.Fit(mipPeakFinderFunc,'QLR0+','',0.4,1.6)
                uncalibMIP,uncalibMIPerr=mipPeakFinderFunc.GetParameter(1),mipPeakFinderFunc.GetParError(1)
                np=uncalibMIPprofiles[cell][iEtaBin].GetN()
                uncalibMIPprofiles[cell][iEtaBin].SetPoint(np,layer,uncalibMIP)
                uncalibMIPprofiles[cell][iEtaBin].SetPointError(np,0,uncalibMIPerr)
                
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
                    dataSent.append( dataBitsToSend(ieta, eH.GetXaxis().GetBinLowEdge(ybin))[0] )
                np=dataVolProfiles[cell][iEtaBin].GetN()
                avgDataVol,stdDataVol=weighted_avg_and_std(dataSent,evCount)
                dataVolProfiles[cell][iEtaBin].SetPoint(np,layer,avgDataVol)
                dataVolProfiles[cell][iEtaBin].SetPointError(np,0,stdDataVol/TMath.Sqrt(totalHits))

                #all done here
                eH.Delete()
                
            #save energy profile canvas
            eprofH.Draw('colz')
            MyPaveText('CMS simulation, <PU>=140 (14 TeV)')
            ptxt=MyPaveText('#bf{[Layer #%d]} cell area=%3.1f mm^{2}'%(layer,baseArea*cell),0.12,0.9,0.9,0.93)
            ptxt.SetTextSize(0.035)
            ptxt.SetFillColor(0)
            ptxt.SetFillStyle(3001)
            eprofGr.Draw('p')
            c.SetRightMargin(0.12)
            c.SetLogz(True)
            c.SetLogy(False)
            c.Modified()
            c.Update()
            c.SaveAs('%s/eprof2D_layer%d_cell%d%s.png'%(output,layer,cell,pfix))

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
            leg.SetHeader('#scale[0.75]{#bf{Pseudo-rapidity}}')
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
                zoomDistsH[ nzoom ].GetYaxis().SetTitle('Events (normalized in range)')
                zoomDistsH[ nzoom ].GetYaxis().SetTitleOffset(1.5)
                totInZoom=zoomDistsH[nzoom].Integral( zoomDistsH[nzoom].GetXaxis().FindBin(0.4),zoomDistsH[nzoom].GetXaxis().FindBin(3) )
                if totInZoom>0 : zoomDistsH[ nzoom ].Scale(1./totInZoom)
                drawOpt='hist same'

            c.cd()
            c.SetLogy()
            c.Modified()
            c.Update()
            c.SaveAs('%s/eprof1D_layer%d_cell%d%s.png'%(output,layer,cell,pfix))

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
            leg.SetHeader('#scale[0.75]{#bf{Pseudo-rapidity}}')
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
            c.SaveAs('%s/occ_thr%d_cell%d%s.png'%(output,10*thr,cell,pfix))
            c.SetLogy(True)
            for eta in occProfiles[thr][cell]: occProfiles[thr][cell][eta].GetYaxis().SetRangeUser(1e-2,100)
            c.Modified()
            c.Update()
            c.SaveAs('%s/occ_thr%d_cell%d%s_log.png'%(output,10*thr,cell,pfix))

    #other summary graphs
    summaryGraphs={'datavol':[dataVolProfiles,    "Average data sent / cell (bit)", [0,3]],
                   'mippeak':[uncalibMIPprofiles, "E / MIP (peak)",                 [0.8,1.4]] 
                   }
    for gr in summaryGraphs:
        for cell in summaryGraphs[gr][0]:
            c.Clear()
            c.SetRightMargin(0.05)
            c.SetLogz(False)
            c.SetLogy(False)
            drawOpt='acp'
            leg=TLegend(0.8,0.72,0.95,0.93)
            leg.SetHeader('#scale[0.75]{#bf{Pseudo-rapidity}}')
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.SetTextFont(42)
            leg.SetTextSize(0.035)
            for eta in summaryGraphs[gr][0][cell]:
                summaryGraphs[gr][0][cell][eta].Draw(drawOpt)
                summaryGraphs[gr][0][cell][eta].GetXaxis().SetTitle("Layer")
                summaryGraphs[gr][0][cell][eta].GetYaxis().SetTitle(summaryGraphs[gr][1])
                summaryGraphs[gr][0][cell][eta].GetYaxis().SetRangeUser(summaryGraphs[gr][2][0],summaryGraphs[gr][2][1])
                drawOpt='cp'
                leg.AddEntry(summaryGraphs[gr][0][cell][eta],summaryGraphs[gr][0][cell][eta].GetTitle(),'cp')

        leg.Draw()
        MyPaveText('CMS simulation, <PU>=140 (14 TeV)')
        MyPaveText('Cell area=%3.1f mm^{2}'%(baseArea*cell),0.12,0.9,0.9,0.93).SetTextSize(0.035)

        c.Modified()
        c.Update()
        c.SaveAs('%s/%s_cell%d%s.png'%(output,gr,cell,pfix))

            
    fIn.Close()

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

    #standard
    drawOccupancyAnalysisSummary(input=opt.input,output=opt.output,mip=opt.mip)

    #using isolation cut (no other deposit >2*MIP in the neighbourhood)
    drawOccupancyAnalysisSummary(input=opt.input,output=opt.output,mip=opt.mip,
                                 addCut='(e_sr1+e_sr2)/16<%f'%opt.mip,
                                 pfix='_iso')

    #adding noise
    drawOccupancyAnalysisSummary(input=opt.input,output=opt.output,mip=opt.mip,
                                 #http://en.wikipedia.org/wiki/Box-Muller_transform
                                 noiseLevel=0.2,
                                 pfix='_noise')
            
if __name__ == "__main__":
    main()

