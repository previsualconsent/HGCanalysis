#!/usr/bin/env python

import ROOT
from UserCode.HGCanalysis.PlotUtils import *
import numpy as np
from array import array
import sys

"""
runs the fits to the calibration and resolution
"""
def showEMcalibrationResults(calibFunc,resFunc):

    resolSummary=[]

    #
    #resolution
    #
    fitFunc=ROOT.TF1('resolmodel',"sqrt([0]*x*x+[1])",0,1)
    #fitFunc=ROOT.TF1('resolmodel',"sqrt([0]*x*x+[1]+[2]*x*x*x*x)",0,1)
    fitFunc.SetParameter(0,0.2);
    fitFunc.SetParLimits(0,0,1);
    fitFunc.SetParameter(1,0);
    fitFunc.SetParLimits(1,0,0.3);

    c=ROOT.TCanvas('cresol','cresol',500,500)
    c.SetTopMargin(0.05)
    leg=ROOT.TLegend(0.15,0.6,0.5,0.93)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    for ir in xrange(0,len(resFunc)):
        if ir==0:
            resFunc[ir].Draw('ap')
            resFunc[ir].GetXaxis().SetTitle("1/#sqrt{E} [1/#sqrt{GeV}]")
            resFunc[ir].GetYaxis().SetTitle("Relative energy resolution") # #sigma_{E}/E
            resFunc[ir].GetXaxis().SetLabelSize(0.04)
            resFunc[ir].GetYaxis().SetLabelSize(0.04)
            resFunc[ir].GetXaxis().SetTitleSize(0.05)
            resFunc[ir].GetYaxis().SetTitleSize(0.05)
            resFunc[ir].GetYaxis().SetTitleOffset(1.3)
            resFunc[ir].GetYaxis().SetRangeUser(0,0.25)
        else:
            resFunc[ir].Draw('p')

        resFunc[ir].Fit(fitFunc,'MER+')
        sigmaStoch=ROOT.TMath.Sqrt(fitFunc.GetParameter(0));
        sigmaStochErr=fitFunc.GetParError(0)/(2*ROOT.TMath.Sqrt(sigmaStoch))
        sigmaConstErr=0
        sigmaConst=ROOT.TMath.Sqrt(fitFunc.GetParameter(1))
        sigmaConstErr=fitFunc.GetParError(1)/(2*ROOT.TMath.Sqrt(sigmaConst))
        #sigmaNoise=ROOT.TMath.Sqrt(fitFunc.GetParameter(2))
        #sigmaNoiseErr=fitFunc.GetParError(2)/(2*ROOT.TMath.Sqrt(sigmaNoise))

        leg.AddEntry(resFunc[ir],
                     #"#splitline{[#bf{#it{%s}}]}{#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.4f#pm%3.4f}"
                     "[#bf{#it{%s}}] #frac{#sigma}{E} #propto #frac{%3.4f #pm %3.4f}{#sqrt{E}} #oplus %3.5f#pm%3.5f"
                     %(resFunc[ir].GetTitle(),sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr),
                     "fp")
        
        resolSummary.append( [sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr] )

    leg.SetTextSize(0.03)
    leg.Draw()
    MyPaveText('CMS simulation')
    c.Modified()
    c.Update()
    #for ext in ['png','pdf','C'] : c.SaveAs(outdir+'/resol.%s'%ext)

    #
    #calibration
    #
    calibModel=ROOT.TF1('calibmodel',"[0]*x+[1]",0,800)
    
    c=ROOT.TCanvas('ccalib','ccalib',500,500)
    c.SetTopMargin(0.05)
    leg=ROOT.TLegend(0.53,0.15,0.9,0.5)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    for ir in xrange(0,len(calibFunc)):
        
        if ir==0:
            calibFunc[ir].Draw('ap')
            calibFunc[ir].GetXaxis().SetTitle("Energy [GeV]")
            calibFunc[ir].GetYaxis().SetTitle("Reconstructed energy")
            calibFunc[ir].GetXaxis().SetLabelSize(0.04)
            calibFunc[ir].GetYaxis().SetLabelSize(0.04)
            calibFunc[ir].GetXaxis().SetTitleSize(0.05)
            calibFunc[ir].GetYaxis().SetTitleSize(0.05)
            calibFunc[ir].GetYaxis().SetTitleOffset(1.3)
        else:
            calibFunc[ir].Draw('p')
        
        calibFunc[ir].Fit(calibModel,'MER+')
        slope=calibModel.GetParameter(0)
        slopeErr=calibModel.GetParError(0)
        offset=calibModel.GetParameter(1)
        offsetErr=calibModel.GetParError(1)

        leg.AddEntry(calibFunc[ir],
                     "#splitline{[#bf{#it{%s}}]}{E #propto (%3.1f#pm%3.1f)#timesE_{rec} +%3.0f#pm%3.0f}"
                     %(calibFunc[ir].GetTitle(),slope,slopeErr,offset,offsetErr),
                     "fp")

        resolSummary[ir].append(slope)
        resolSummary[ir].append(slopeErr)
        resolSummary[ir].append(offset)
        resolSummary[ir].append(offsetErr)
        
    leg.Draw()
    MyPaveText('CMS simulation')
    c.Modified()
    c.Update()
    #for ext in ['png','pdf','C'] : c.SaveAs(outdir+'/calib.%s'%ext)

    return resolSummary

customROOTstyle()
ROOT.gROOT.SetBatch(False)

ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

"""
"""
def runResolutionStudy(vars=[["recEn","#Sigma w_{i}E_{i}"],["recEnDR1","#Sigma_{#Delta R<0.1}w_{i}E_{i}"]],
                       wgt=[1.,3.3/1.6,5.6/1.6],
                       mipEn=55.1,
                       url='cmssw/SingleElectron_v15.root'):

    #prepare the workspace
    ws=ROOT.RooWorkspace("w")
    dsVars=ROOT.RooArgSet( ws.factory('pt[0,0,9999999999]'), ws.factory('en[0,0,9999999999]') )
    varNames=[]
    for var,varTitle in vars:
        dsVars.add( ws.factory('%s[0,0,9999999999]'%var) )
        varNames.append(var)
    getattr(ws,'import')( ROOT.RooDataSet('data','data',dsVars) )

    binGenEn = [10,20,30,40,50,75,100,150,200,250,400] #,400,600,1000]
    hen=ROOT.TH1F('hgen',';Generated energy [GeV];Events',len(binGenEn)-1,array('f',binGenEn))
    hen.SetDirectory(0)

    #read all to a RooDataSet
    fin=ROOT.TFile.Open('cmssw/SingleElectron_v15.root')
    HGC=fin.Get('HGC')
    for entry in xrange(0,HGC.GetEntriesFast()+1):

        HGC.GetEntry(entry)
        
        if ROOT.TMath.Abs(HGC.genEta)<1.6 or ROOT.TMath.Abs(HGC.genEta)>2.9 : continue
        #if ROOT.TMath.Abs(HGC.genEta)<2 or ROOT.TMath.Abs(HGC.genEta)>2.2 : continue
        #if ROOT.TMath.Abs(HGC.genEta)<2.5 or ROOT.TMath.Abs(HGC.genEta)>2.6 : continue
        #if ROOT.TMath.Abs(HGC.genEta)<2.8 or ROOT.TMath.Abs(HGC.genEta)>2.9 : continue
            
        ws.var('pt').setVal(HGC.genPt)
        ws.var('en').setVal(HGC.genPt*ROOT.TMath.CosH(HGC.genEta))
        hen.Fill( ws.var('en').getVal() )
        newEntry=ROOT.RooArgSet( ws.var('pt'), ws.var('en') )

        #energy estimators
        enCorr=(ROOT.TMath.TanH(HGC.genEta)*1e6/mipEn)/ws.var('en').getVal()
        if 'recEn' in varNames:
            ws.var('recEn').setVal(     (wgt[0]*HGC.ec1_en     + wgt[1]*HGC.ec2_en     + wgt[2]*HGC.ec3_en)*enCorr )
            newEntry.add( ws.var('recEn') )
        if 'recEnDR1' in varNames:
            ws.var('recEnDR1').setVal(  (wgt[0]*HGC.ec1_endr1  + wgt[1]*HGC.ec2_endr1  + wgt[2]*HGC.ec3_endr1 )*enCorr )
            newEntry.add( ws.var('recEnDR1') )
        if 'recEnDR25' in varNames:
            ws.var('recEnDR25').setVal( (wgt[0]*HGC.ec1_endr25 + wgt[1]*HGC.ec2_endr25 + wgt[3]*HGC.ec3_endr25)*enCorr )
            newEntry.add( ws.var('recEnDR25') )
        if 'recEnDR5' in varNames:
            ws.var('recEnDR5').setVal(  (wgt5[0]*HGC.ec1_endr5  + wgt[1]*HGC.ec2_endr5  + wgt[3]*HGC.ec3_endr5 )*enCorr )
            newEntry.add( ws.var('recEnDR5') )
        ws.data('data').add( newEntry )
    fin.Close()

     
    #prepare to fit the energy slices
    calibFunc=[]
    resFunc=[]
    c=ROOT.TCanvas('c','c',1400,800)
    for var,varTitle in vars:
        nv=len(calibFunc)
        calibFunc.append(ROOT.TGraphErrors())
        calibFunc[nv].SetMarkerStyle(20+nv)
        calibFunc[nv].SetTitle(varTitle)
        calibFunc[nv].SetFillStyle(0)
        calibFunc[nv].SetName('calib_'+var)    
        resFunc.append(calibFunc[nv].Clone('resol_'+var))

        #run resolution fits in different energy ranges
        c.Clear()
        c.Divide(hen.GetXaxis().GetNbins()/3+1,3)
        ipad=0
        for ibin in xrange(1,hen.GetXaxis().GetNbins()+1):
            ipad=ipad+1
            postfix='fit%d'%ipad
            enmin=hen.GetXaxis().GetBinLowEdge(ibin)
            enmax=hen.GetXaxis().GetBinUpEdge(ibin)
            redData = ws.data('data').reduce("en>=%f && en<=%f"%(enmin,enmax))
            if redData.sumEntries()<10 : continue
            redDataMean=redData.mean(ws.var(var))
            redDataSigma=redData.sigma(ws.var(var))
            ws.var(var).setRange(postfix,redDataMean-6*redDataSigma,redDataMean+6*redDataSigma)
            ws.var(var).setRange('fit'+postfix,redDataMean-1.5*redDataSigma,redDataMean+1.5*redDataSigma)
            ws.factory('Gaussian::g_%s(%s,mean_%s[%f,%f,%f],sigma_%s[%f,%f,%f])'%
                       (postfix,var,
                        postfix,redDataMean,redDataMean-5*redDataSigma,redDataMean+5*redDataSigma,
                           postfix,redDataSigma,redDataSigma/2,redDataSigma*2)
                )

            #fit
            fres=ws.pdf('g_%s'%postfix).fitTo( redData,ROOT.RooFit.Range('fit'+postfix),ROOT.RooFit.Save(True))
            meanFit=ws.var('mean_%s'%postfix).getVal()
            meanFitError=ws.var('mean_%s'%postfix).getError()
            sigmaFit=ws.var('sigma_%s'%postfix).getVal()
            sigmaFitError=ws.var('sigma_%s'%postfix).getError()

            np=calibFunc[nv].GetN()
            meanGenEn=0.5*(enmax+enmin)
            meanGenEnErr=0.5*(enmax-enmin)
            calibFunc[nv].SetPoint(np,meanGenEn,meanFit*meanGenEn)
            calibFunc[nv].SetPointError(np,meanGenEnErr,meanFitError*meanGenEn)
            resFunc[nv].SetPoint(np,1/ROOT.TMath.Sqrt(meanGenEn),sigmaFit/meanFit)
            resFunc[nv].SetPointError(np,meanGenEnErr/(2*ROOT.TMath.Power(meanGenEn,1.5)),sigmaFitError/meanFit)
            
            #show the result
            p=c.cd(ipad)
            frame=ws.var(var).frame(ROOT.RooFit.Range(postfix))
            redData.plotOn(frame)
            ws.pdf('g_%s'%postfix).plotOn(frame,ROOT.RooFit.Range(postfix))
            frame.Draw()
            frame.GetXaxis().SetNdivisions(5)
            frame.GetXaxis().SetTitle(varTitle + ' / Energy')
            #frame.GetYaxis().SetRangeUser(0,1600)
            frame.GetYaxis().SetTitle('Events')
            frame.GetYaxis().SetTitleOffset(1.4)
            pt=MyPaveText('[%d<=Energy/GeV<%d]\\<E>=%3.1f RMS=%3.1f\\#mu=%3.1f #sigma=%3.1f'%(enmin,enmax,redDataMean,redDataSigma,meanFit,sigmaFit),
                          0.18,0.9,0.5,0.6)
            pt.SetTextFont(42)
            pt.SetTextSize(0.06)
            if ipad==1:
                pt=MyPaveText('CMS simulation')
                pt.SetTextSize(0.08)

    return showEMcalibrationResults(calibFunc=calibFunc,resFunc=resFunc)

"""
steer 
"""
def main(argv=None):

    customROOTstyle()
    ROOT.gROOT.SetBatch(False)

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    w2Scan=[0.5,0.75,1.0,1.5,2.0,2.5,3.0,3.5]
    w3Scan=[0.5,0.75,1.0,1.5,2.0,2.5,3.0,3.5]
    hstoch=ROOT.TH2F('stochscan',';w_{2}/w_{1};w_{3}/w_{1};Value',len(w2Scan),0,len(w2Scan),len(w3Scan),0,len(w3Scan))
    for xbin in xrange(1,hstoch.GetXaxis().GetNbins()+1): hstoch.GetXaxis().SetBinLabel(xbin,'%3.2f'%w2Scan[xbin-1])
    for ybin in xrange(1,hstoch.GetYaxis().GetNbins()+1): hstoch.GetYaxis().SetBinLabel(ybin,'%3.2f'%w3Scan[ybin-1])
    hstoch.SetTitle('(#sigma/E)_{stoch}')
    hconst=hstoch.Clone('constscan')
    hconst.SetTitle('(#sigma/E)_{const}')
    hslope=hstoch.Clone('slopescan')
    hslope.SetTitle('Calibration slope')
    hoffset=hstoch.Clone('offsetscan')
    hoffset.SetTitle('Calibration offset')
    hstoch.SetDirectory(0)
    hconst.SetDirectory(0)
    hslope.SetDirectory(0)
    hoffset.SetDirectory(0)
    xbin=0
    for w2 in w2Scan:
        xbin=xbin+1
        ybin=0
        for w3 in w3Scan:
            ybin=ybin+1
            res=runResolutionStudy(vars=[['recEn','#Sigma w_{i}E_{i}']],
                                   wgt=[1.0,w2,w3],
                                   mipEn=55.1,
                                   url='cmssw/SingleElectron_v15.root')
            hstoch.SetBinContent(xbin,ybin,res[0][0])
            hstoch.SetBinError(xbin,ybin,res[0][1])
            hconst.SetBinContent(xbin,ybin,res[0][2])
            hconst.SetBinError(xbin,ybin,res[0][3])
            hslope.SetBinContent(xbin,ybin,res[0][4])
            hslope.SetBinError(xbin,ybin,res[0][5])
            hoffset.SetBinContent(xbin,ybin,res[0][6])
            hoffset.SetBinError(xbin,ybin,res[0][7])

    #show results
    c=ROOT.TCanvas('cscan','cscan',600,600)
    c.SetRightMargin(0.15)
    for h in [hstoch,hconst,hslope,hoffset]:
        c.Clear()
        h.Draw('colz')
        h.GetXaxis().SetLabelSize(0.035)
        h.GetYaxis().SetLabelSize(0.035)
        h.GetXaxis().SetTitleSize(0.04)
        h.GetYaxis().SetTitleSize(0.04)
        MyPaveText('CMS simulation')
        pt=MyPaveText(h.GetTitle(),0.18,0.9,0.5,0.6)
        pt.SetTextFont(42)
        pt.SetTextSize(0.06)
        c.Modified()
        c.Update()
        raw_input()
        c.SaveAs('%s.png'%h.GetName())

    
if __name__ == "__main__":
    sys.exit(main())
