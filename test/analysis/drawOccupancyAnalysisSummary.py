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
Return the weighted average and standard deviation.values, weights -- Numpy ndarrays with the same shape.
"""
def weighted_avg_and_std(values, weights):
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return average, TMath.Sqrt(variance)

"""
"""
def drawOccupancyAnalysisSummary(input,output,mip,addCut='',noiseLevel=0,pfix='',caption='') :

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
    baseArea=1
    for ientry in xrange(0,occT.GetEntriesFast()-1):
        occT.GetEntry(ientry)
        if occT.cell!=1 : continue 
        baseArea=occT.area
        break

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

    #base canvas
    c=TCanvas('c','c',500,500)
    
    cells=[2]
    #thresholds=[0.4]
    #cells=[1,2,3]
    thresholds=[0.4,1.0,5.0,10.,25.0]
    etaBins=[[1,3],[4,5],[6,7],[8,9],[10,11],[12,13]]#,[14,14],[15,15]]

    #generate code to compress data
    #from ROOT import testCompressionAlgos, getTriggerBits, getReadoutBits
    #roespectrumH=TH1F('roespectrumH',';(E/MIP) / coth(#eta);Hits (all layers)',255,0.4,102.4)
    #roespectrumH.Sumw2()
    #roespectrumH.SetLineWidth(2)
    #trigespectrumH=TH1F('trigespectrumH',';(E/MIP) / coth(#eta);Hits (all layers)',255,10,1285)
    #trigespectrumH.Sumw2()
    #trigespectrumH.SetLineWidth(2)
    #for etaRange in [[1.5,2.0],[2.0,2.5],[2.5,3.0]]:
    #    for cell in cells:
    #        roespectrumH.Reset('ICE')
    #        occT.Draw('(e_sr0/%f)*TMath::TanH(eta) >> roespectrumH'%mip,'cell==%d && TMath::Abs(eta)<%f && TMath::Abs(eta)>=%f'%(cell,etaRange[1],etaRange[0]),'goff')
    #        trigespectrumH.Reset('ICE')
    #        occT.Draw('(e_sr0/%f)*TMath::TanH(eta) >> trigespectrumH'%mip,'cell==%d && TMath::Abs(eta)<%f && TMath::Abs(eta)>=%f'%(cell,etaRange[1],etaRange[0]),'goff')
    #        res=testCompressionAlgos(roespectrumH,trigespectrumH,0.5*(etaRange[1]+etaRange[0]))

    #        report='Difference wrt to baseline (%)'
    #        for r in res: report = report +'\\%s %3.1f'%(r.first,100*(r.second-1))
        
            #save energy profile canvas
    #        c.Clear()
    #        roespectrumH.Draw('hist')
    #        MyPaveText(caption)
    #        ptxt=MyPaveText('#bf{%3.1f<|#eta|<%3.1f}\\Cell area=%3.1f mm^{2}\\%s'%(etaRange[0],etaRange[1],baseArea*cell,report),0.4,0.65,0.95,0.93)
    #        ptxt.SetTextSize(0.035)
    #        ptxt.SetFillColor(0)
    #        ptxt.SetFillStyle(3001)
    #        c.SetRightMargin(0.12)
    #        c.SetLogy(True)
    #        c.Update()
    #        c.SaveAs('%s/espectrum_eta%dto%d_cell%d%s.png'%(output,etaRange[0]*10,etaRange[1]*10,cell,pfix))

    #prepare plots
    eprofH=TH2F('eprofH',';Pseudo-rapidity;(E/MIP) / coth(#eta);Hits',15,1.5,2.9,130,0,26)
    eprofH.GetZaxis().SetTitleOffset(-0.3)
    eprofH.GetZaxis().SetLabelSize(0.035)
    eprofH.Sumw2()
    eprofGr=TGraphErrors()
    eprofGr.SetMarkerStyle(20)
    eprofGr.SetMarkerColor(38)
    eprofGr.SetLineColor(38)
    eprofGr.SetLineWidth(2)
    occProfiles={}
    dataVolH=TH2F('datavolH',';Pseudo-rapidity;#bits sent/bunch crossing;Hits',15,1.5,2.9,32,0,32)
    dataVolH.GetZaxis().SetTitleOffset(-0.3)
    dataVolH.GetZaxis().SetLabelSize(0.035)
    dataVolH.Sumw2()
    dataprofGr=TGraphErrors()
    dataprofGr.SetMarkerStyle(20)
    dataprofGr.SetMarkerColor(38)
    dataprofGr.SetLineColor(38)
    dataprofGr.SetLineWidth(2)
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
       
    #iterate over layers
    for layer in layers:

        if layer%3==2 or layer%3==0 : continue

        sys.stdout.write( '\r Layer [%d/%d]'%(layer/3+1,len(layers)/3) )
        sys.stdout.flush()

        for cell in cells:

            #create an entry list with these cuts
            finalCut='(layer==%d || layer==%d) && cell==%d'%(layer,layer+1,cell)
            if len(addCut) : finalCut=finalCut + ' && ' + addCut            
            occT.Draw('>>elist', finalCut,'entrylist')
            
            elist = ROOT.gDirectory.Get("elist")

            #fill the base histograms 
            dataVolH.Reset('ICE')
            eprofH.Reset('ICE')
            edeps={}
            for ientry in xrange(0,elist.GetN()):
                tentry=elist.GetEntry(ientry) 
                occT.GetEntry( tentry )

                corrEn=(occT.e_sr0/mip)*TMath.TanH(occT.eta)
                if occ.layer%3==1: edeps[occT.eta]=corrEn
                else : edeps[occT.eta]=corrEn
            for eta in edeps:
                corrEn=edeps[eta]
                if noiseLevel>0 :
                    corrEn = corrEn + gRandom.Gaus(0,noiseLevel)
                    if corrEn<0: corrEn=0
                dataVolH.Fill(TMath.Abs(eta),getReadoutBits(corrEn,eta))
                eprofH.Fill(TMath.Abs(eta),corrEn)

            #project in eta slices
            eprofGr.Set(0)
            dataprofGr.Set(0)
            eDistsH=[]
            dataDistsH=[]
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
                meanEn,meanEnErr=eH.GetMean(),eH.GetRMS()
                np=eprofGr.GetN()
                eprofGr.SetPoint(np,ieta,meanEn)
                eprofGr.SetPointError(np,0.5*deta,meanEnErr)
                    
                #fit near the MIP peak
                mipPeakFinderFunc.SetParameter(1,1)
                eH.Fit(mipPeakFinderFunc,'QLR0+','',0.4,1.6)
                uncalibMIP,uncalibMIPerr=mipPeakFinderFunc.GetParameter(1),mipPeakFinderFunc.GetParError(1)
                np=uncalibMIPprofiles[cell][iEtaBin].GetN()
                uncalibMIPprofiles[cell][iEtaBin].SetPoint(np,layer/3+1,uncalibMIP)
                uncalibMIPprofiles[cell][iEtaBin].SetPointError(np,0,uncalibMIPerr)
                
                #compute occupancies
                for thr in thresholds:
                    ybin=eH.GetXaxis().FindBin(thr)
                    cellOccErr=ROOT.Double(0)
                    cellOcc=eH.IntegralAndError(ybin,eH.GetXaxis().GetNbins(),cellOccErr)
                    np=occProfiles[thr][cell][iEtaBin].GetN()
                    occProfiles[thr][cell][iEtaBin].SetPoint(np,layer/3+1,cellOcc*100/totalHits)
                    occProfiles[thr][cell][iEtaBin].SetPointError(np,0,cellOccErr*100/totalHits)
                    
                #data volume
                dvH=dataVolH.ProjectionY('dvH',etaBins[iEtaBin][0],etaBins[iEtaBin][1])
                dvH.SetTitle('%3.1f-%3.1f'%(etaMin,etaMax))
                fixExtremities(dvH)
                dataDistsH.append( dvH.Clone('dvprof_%d'%iEtaBin) )
                meanData,meanDataErr=dvH.GetMean(),dvH.GetRMS()
                np=dataprofGr.GetN()
                dataprofGr.SetPoint(np,ieta,meanData)
                dataprofGr.SetPointError(np,0.5*deta,meanDataErr)
                np=dataVolProfiles[cell][iEtaBin].GetN()
                dataVolProfiles[cell][iEtaBin].SetPoint(np,layer/3+1,meanData)
                dataVolProfiles[cell][iEtaBin].SetPointError(np,0,dvH.GetMeanError())

                #all done here
                eH.Delete()
                dvH.Delete()
                
            #save profile canvas
            profilesToSave={'eprof2D':[eprofH,eprofGr],
                            'dvprof2D':[dataVolH,dataprofGr]}
            for prof in profilesToSave:
                c.Clear()
                profilesToSave[prof][0].Draw('colz')
                profilesToSave[prof][1].Draw('p')
                MyPaveText(caption)
                ptxt=MyPaveText('#bf{[Layer #%d]} cell area=%3.1f mm^{2}'%(layer/3+1,baseArea*cell),0.12,0.9,0.9,0.93)
                ptxt.SetTextSize(0.035)
                ptxt.SetFillColor(0)
                ptxt.SetFillStyle(3001)
                c.SetRightMargin(0.12)
                c.SetLogz(True)
                c.SetLogy(False)
                c.Modified()
                c.Update()
                c.SaveAs('%s/%s_layer%d_cell%d%s.png'%(output,prof,layer/3+1,cell,pfix))
                c.SaveAs('%s/%s_layer%d_cell%d%s.C'%(output,prof,layer/3+1,cell,pfix))

            #distributions at a given #eta
            slicesToSave={'eprof1D':[eDistsH,[0.4,3]],
                          'dv1D':   [dataDistsH,[0,16]]}
            for dist in slicesToSave:
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
                for ih in xrange(len(slicesToSave[dist][0])-1,0,-2):
                    slicesToSave[dist][0][ih].Draw(drawOpt)
                    drawOpt='hist same'
                    slicesToSave[dist][0][ih].GetYaxis().SetTitle('Hits')
                    ctr=(len(slicesToSave[dist][0])-ih)/2
                    slicesToSave[dist][0][ih].SetLineColor(1+ctr/2)
                    slicesToSave[dist][0][ih].SetLineStyle(1+ctr%2)
                    slicesToSave[dist][0][ih].SetLineWidth(2)
                    leg.AddEntry(slicesToSave[dist][0][ih],slicesToSave[dist][0][ih].GetTitle(),'l')
                leg.Draw()
                MyPaveText(caption)
                ptxt=MyPaveText('#bf{[Layer #%d]} cell area=%3.1f mm^{2}'%(layer/3+1,baseArea*cell),0.12,0.9,0.9,0.93)
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
                for ih in xrange(len(slicesToSave[dist][0])-1,0,-2):
                    nzoom=len(zoomDistsH)
                    zoomDistsH.append( slicesToSave[dist][0][ih].Clone('%s_zoom'%slicesToSave[dist][0][ih].GetName()) )
                    zoomDistsH[ nzoom ].Draw(drawOpt)
                    zoomMin=slicesToSave[dist][1][0]
                    zoomMax=slicesToSave[dist][1][1]
                    zoomDistsH[ nzoom ].GetXaxis().SetRangeUser(zoomMin,zoomMax)
                    zoomDistsH[ nzoom ].GetYaxis().SetTitle('Events (normalized in range)')
                    zoomDistsH[ nzoom ].GetYaxis().SetTitleOffset(1.5)
                    totInZoom=zoomDistsH[nzoom].Integral( zoomDistsH[nzoom].GetXaxis().FindBin(zoomMin),zoomDistsH[nzoom].GetXaxis().FindBin(zoomMax) )
                    if totInZoom>0 : zoomDistsH[ nzoom ].Scale(1./totInZoom)
                    drawOpt='hist same'
                    
                c.cd()
                c.SetLogy()
                c.Modified()
                c.Update()
                c.SaveAs('%s/%s_layer%d_cell%d%s.png'%(output,dist,layer/3+1,cell,pfix))
                c.SaveAs('%s/%s_layer%d_cell%d%s.C'%(output,dist,layer/3+1,cell,pfix))

    #occupancy summaries
    for thr in occProfiles:
        for cell in occProfiles[thr]:

            c.Clear()
            c.SetRightMargin(0.05)
            c.SetLogz(False)
            c.SetLogy(False)
            drawOpt='acp'
            leg=TLegend(0.75,0.72,0.95,0.93)
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
            MyPaveText(caption)
            MyPaveText('Cell area=%3.1f mm^{2}, Threshold=%3.1f/MIP'%(baseArea*cell,thr),0.12,0.9,0.9,0.93).SetTextSize(0.035)

            c.Modified()
            c.Update()
            c.SaveAs('%s/occ_thr%d_cell%d%s.png'%(output,10*thr,cell,pfix))
            c.SaveAs('%s/occ_thr%d_cell%d%s.C'%(output,10*thr,cell,pfix))
            c.SetLogy(True)
            for eta in occProfiles[thr][cell]: occProfiles[thr][cell][eta].GetYaxis().SetRangeUser(1e-2,100)
            c.Modified()
            c.Update()
            c.SaveAs('%s/occ_thr%d_cell%d%s_log.png'%(output,10*thr,cell,pfix))
            c.SaveAs('%s/occ_thr%d_cell%d%s_log.C'%(output,10*thr,cell,pfix))

    #other summary graphs
    summaryGraphs={'datavol':[dataVolProfiles,    "Average data sent / cell (bit)", [0,10]],
                   'mippeak':[uncalibMIPprofiles, "E / MIP (peak)",                 [0.5,2.0]] 
                   }
    for gr in summaryGraphs:
        for cell in summaryGraphs[gr][0]:
            c.Clear()
            c.SetRightMargin(0.05)
            c.SetLogz(False)
            c.SetLogy(False)
            drawOpt='ap'
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
                drawOpt='p'
                leg.AddEntry(summaryGraphs[gr][0][cell][eta],summaryGraphs[gr][0][cell][eta].GetTitle(),'cp')

            leg.Draw()
            MyPaveText(caption)
            MyPaveText('Cell area=%3.1f mm^{2}'%(baseArea*cell),0.12,0.9,0.9,0.93).SetTextSize(0.035)

            c.Modified()
            c.Update()
            c.SaveAs('%s/%s_cell%d%s.png'%(output,gr,cell,pfix))
            c.SaveAs('%s/%s_cell%d%s.C'%(output,gr,cell,pfix))
            
    fIn.Close()

"""
checks the input arguments and steers the analysis
"""
def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in'       ,    dest='input'          , help='Input file',       default=None)
    parser.add_option('-o', '--out'      ,    dest='output'         , help='Output dir',       default=None)
    parser.add_option('-m', '--mip'      ,    dest='mip'            , help='MIP energy [keV]', default=55.1, type=float)
    parser.add_option('-c', '--caption'  ,    dest='caption'        , help='Global caption',   default='CMS simulation, <PU>=140 (14 TeV)')
    (opt, args) = parser.parse_args()

    if opt.input is None or opt.output is None:
        parser.print_help()
        sys.exit(1)

    #plot formatting
    customROOTstyle()
    gROOT.SetBatch(True)
    #gROOT.SetBatch(False)
    gStyle.SetPalette(53)

    #standard
    drawOccupancyAnalysisSummary(input=opt.input,output=opt.output,mip=opt.mip,caption=opt.caption)

    #using isolation cut (no other deposit >2*MIP in the neighbourhood)
    drawOccupancyAnalysisSummary(input=opt.input,output=opt.output,mip=opt.mip,
                                 addCut='(e_sr1+e_sr2)/16<%f'%opt.mip,
                                 pfix='_iso',
                                 caption=opt.caption)
    
    #adding noise
    drawOccupancyAnalysisSummary(input=opt.input,output=opt.output,mip=opt.mip,
                                 #http://en.wikipedia.org/wiki/Box-Muller_transform
                                 noiseLevel=1./5.,
                                 pfix='_noise',
                                 caption=opt.caption)
            
if __name__ == "__main__":
    main()


