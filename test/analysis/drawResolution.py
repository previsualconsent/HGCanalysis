#!/usr/bin/env python

import ROOT
from PlotUtils import *
import numpy as numpy
from array import array
import io,os,sys
import json
    
"""
Converts all to a workspace and returns optimized weights
"""
def prepareWorkspace(url,integSD,integRanges,mipEn,genM,outUrl):
    
    #prepare the workspace
    ws=ROOT.RooWorkspace("w")
    dsVars=ROOT.RooArgSet( ws.factory('pt[0,0,9999999999]'),  ws.factory('eta[1.5,1.45,3.1]'), ws.factory('en[0,0,9999999999]') )
    for ireg in xrange(0,len(integRanges)): dsVars.add( ws.factory('edep%d[0,0,99999999]'%ireg) )
    getattr(ws,'import')( ROOT.RooDataSet('data','data',dsVars) )

    #optimization with linear regression
    optimVec    = numpy.zeros( len(integRanges)+1 )
    optimMatrix = numpy.zeros( (len(integRanges)+1, len(integRanges)+1 ) )

    #read all to a RooDataSet
    fin=ROOT.TFile.Open(url)
    HGC=fin.Get('HGC')
    for entry in xrange(0,HGC.GetEntriesFast()+1):

        HGC.GetEntry(entry)

        genEta=ROOT.TMath.Abs(HGC.genEta)
        if genEta<1.5 or genEta>2.9 : continue
        ws.var('eta').setVal(genEta)

        genPt=HGC.genPt
        genPl=genPt*ROOT.TMath.SinH(genEta)
        ws.var('pt').setVal(genPt)
        genP=genPt*ROOT.TMath.CosH(genEta)
        genEn=ROOT.TMath.Sqrt(genP*genP+genM*genM)
        ws.var('en').setVal(genEn)
        genY=0.5*ROOT.TMath.Log((genEn+genPl)/(genEn-genPl))

        newEntry=ROOT.RooArgSet( ws.var('pt'),  ws.var('eta'), ws.var('en') )

        #get the relevant energy deposits and add new row
        edeps=[]
        for ireg in xrange(0,len(integRanges)):
            sdPrefix='s%dp'%integSD[ireg]
            if HGC.genEta<0 : sdPrefix='s%dm'%integSD[ireg]
            totalEnInIntegRegion=0
            for ilayer in xrange(integRanges[ireg][0],integRanges[ireg][1]+1):
                enCorr=ROOT.TMath.Abs(ROOT.TMath.TanH(genEta)*1e6/mipEn[ireg])
                totalEnInIntegRegion=totalEnInIntegRegion+(getattr(HGC,sdPrefix+'%d_en'%ilayer))*enCorr
            edeps.append(totalEnInIntegRegion)
            ws.var('edep%d'%ireg).setVal(totalEnInIntegRegion)
            newEntry.add(ws.var('edep%d'%(ireg)))

        #a mip veto
        #if ws.var('edep0').getVal()>1 : continue
            
        ws.data('data').add( newEntry )

        #fill the optmization matrix and vector
        for ie in xrange(0,len(edeps)):
            optimVec[ie]=optimVec[ie]+edeps[ie]*genEn
            for je in xrange(0,len(edeps)):
                optimMatrix[ie][je]=optimMatrix[ie][je]+edeps[ie]*edeps[je]
            optimMatrix[len(edeps)][ie]=optimMatrix[len(edeps)][ie]+edeps[ie]
            optimMatrix[ie][len(edeps)]=optimMatrix[ie][len(edeps)]+edeps[ie]
        optimMatrix[len(edeps)][len(edeps)]=optimMatrix[len(edeps)][len(edeps)]+1
        optimVec[len(edeps)]=optimVec[len(edeps)]+genEn
        
    fin.Close()

    #all done, write to file
    ws.writeToFile(outUrl,True)

    #finalize weight optimization
    optimWeights=numpy.linalg.solve(optimMatrix,optimVec)

    #all done here
    return optimWeights

    

"""
runs the fits to the calibration and resolution
"""
def showEMcalibrationResults(calibFunc,resFunc,outDir):

    resolSummary=[]

    #
    #resolution
    #
    fitFunc=ROOT.TF1('resolmodel',"sqrt([0]*x*x+[1])",0,1)
    #fitFunc=ROOT.TF1('resolmodel',"sqrt([0]*x*x+[1]+[2]*x*x*x*x)",0,1)
    fitFunc.SetParameter(0,0.2);
    fitFunc.SetParLimits(0,0,2);
    fitFunc.SetParameter(1,0);
    fitFunc.SetParLimits(1,0,1.0);
    fitFunc.SetLineWidth(1)

    c=ROOT.TCanvas('cresol','cresol',1200,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.5)
    leg=ROOT.TLegend(0.52,0.1,0.99,0.2+0.085*len(calibFunc)/3)
    leg.SetNColumns(3)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.025)
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
            maxY=resFunc[ir].GetMaximum()
            resFunc[ir].GetYaxis().SetRangeUser(0,1.0)
            #resFunc[ir].GetYaxis().SetRangeUser(0,0.1)
        else:
            resFunc[ir].Draw('p')

        resFunc[ir].Fit(fitFunc,'MER+')
        resFunc[ir].GetFunction(fitFunc.GetName()).SetLineColor(resFunc[ir].GetLineColor())
        sigmaStochErr=0
        sigmaStoch=ROOT.TMath.Sqrt(fitFunc.GetParameter(0));
        if sigmaStoch>0 : sigmaStochErr=fitFunc.GetParError(0)/(2*ROOT.TMath.Sqrt(sigmaStoch))
        sigmaConstErr=0
        sigmaConst=ROOT.TMath.Sqrt(fitFunc.GetParameter(1))
        if sigmaConst>0 : sigmaConstErr=fitFunc.GetParError(1)/(2*ROOT.TMath.Sqrt(sigmaConst))
        #sigmaNoise=ROOT.TMath.Sqrt(fitFunc.GetParameter(2))
        #sigmaNoiseErr=fitFunc.GetParError(2)/(2*ROOT.TMath.Sqrt(sigmaNoise))

        leg.AddEntry(resFunc[ir],
                     #"#splitline{[#bf{#it{%s}}]}{#frac{#sigma}{E} #propto #frac{%3.3f #pm %3.3f}{#sqrt{E}} #oplus %3.4f#pm%3.4f}"
                     #"[#bf{#it{%s}}] #frac{#sigma}{E} #propto #frac{%3.4f #pm %3.4f}{#sqrt{E}} #oplus %3.5f#pm%3.5f"
                     "#splitline{[#scale[0.8]{#bf{#it{%s}}}]}{#frac{#sigma}{E} #propto #frac{%3.4f}{#sqrt{E}} #oplus %3.5f}"
                     %(resFunc[ir].GetTitle(),sigmaStoch,sigmaConst),
                     "fp")
        
        resolSummary.append( [sigmaStoch,sigmaStochErr,sigmaConst,sigmaConstErr] )

    leg.SetTextSize(0.03)
    leg.Draw()
    MyPaveText('CMS simulation')
    c.Modified()
    c.Update()
    for ext in ['png','pdf','C'] : c.SaveAs('%s/resolution.%s'%(outDir,ext))

    
    #
    #calibration
    #
    calibModel=ROOT.TF1('calibmodel',"[0]*x+[1]",0,800)
    
    c=ROOT.TCanvas('ccalib','ccalib',1200,500)
    c.SetTopMargin(0.05)
    c.SetRightMargin(0.5)
    leg=ROOT.TLegend(0.52,0.1,0.99,0.2+0.085*len(calibFunc)/3)
    leg.SetNColumns(3)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.025)
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
        calibFunc[ir].GetFunction(calibModel.GetName()).SetLineColor(calibFunc[ir].GetLineColor())
        slope=calibModel.GetParameter(0)
        slopeErr=calibModel.GetParError(0)
        offset=calibModel.GetParameter(1)
        offsetErr=calibModel.GetParError(1)

        leg.AddEntry(calibFunc[ir],
                     #"#splitline{[#bf{#it{%s}}]}{E #propto (%3.1f#pm%3.1f)#timesE_{rec} +%3.0f#pm%3.0f}"
                     "#splitline{[#scale[0.8]{#bf{#it{%s}}}]}{#hat{E} = %3.1f#timesE_{beam} + %3.0f}"
                     %(calibFunc[ir].GetTitle(),slope,offset),
                     "fp")

        resolSummary[ir].append(slope)
        resolSummary[ir].append(slopeErr)
        resolSummary[ir].append(offset)
        resolSummary[ir].append(offsetErr)
        
    leg.Draw()
    MyPaveText('CMS simulation')
    c.Modified()
    c.Update()
    for ext in ['png','pdf','C'] : c.SaveAs('%s/calibration.%s'%(outDir,ext))

    return resolSummary


"""
"""
def runResolutionStudy(url,genM):

    mipEn      = [55.1,        55.1,        55.1,        55.1,        85.0,     1498.4]   
    integSD    = [0,           0,           0,           0,           1,        2]
    integRanges= [[1,1],       [2,11],      [12,21],     [22,31],     [1,12],   [1,10]]
    #defWeights = [0.209/0.494, 0.494/0.494, 0.797/0.494, 1.207/0.494, 0.,       0.,      0.]
    defWeights = [0.028,       0.065,        0.105,       0.160,      1.0,       1.667,   0.]




    #read particle gun configuration
    pGunF=ROOT.TFile.Open('particle_gun_pdf.root')
    pGunPt=pGunF.Get("pt")
    pGunPtRanges=[]
    for np in xrange(0,pGunPt.GetN()):
        x,y=ROOT.Double(0),ROOT.Double(0)
        pGunPt.GetPoint(np,x,y)
        if y>0: 
            pGunPtRanges.append([x])
            nrang=len(pGunPtRanges)
            if nrang==1: continue
            dPt=0.25*(pGunPtRanges[nrang-1][0]-pGunPtRanges[nrang-2][0])
            pGunPtRanges[nrang-1].append(dPt)
            if nrang==2: pGunPtRanges[0].append(dPt)
    pGunY=pGunF.Get("rapidity")
    pGunYRanges=[]
    for np in xrange(0,pGunY.GetN()):
        x,y=ROOT.Double(0),ROOT.Double(0)
        pGunY.GetPoint(np,x,y)
        if y>0: 
            pGunYRanges.append([x])
            nrang=len(pGunYRanges)
            if nrang==1: continue
            dY=0.25*(pGunYRanges[nrang-1][0]-pGunYRanges[nrang-2][0])
            pGunYRanges[nrang-1].append(dY)
            if nrang==2: pGunYRanges[0].append(dY)
    pGunF.Close()

    #prepare output
    outDir=url.replace('.root','')
    os.system('mkdir -p '+outDir)

    #get workspace
    wsOutUrl=url.replace('.root','_ws.root')    
    optimWeights=prepareWorkspace(url=url,integSD=integSD,integRanges=integRanges,mipEn=mipEn,genM=genM,outUrl=url.replace('.root','_ws.root'))
    wsOutF=ROOT.TFile.Open(wsOutUrl)
    ws=wsOutF.Get('w')
    wsOutF.Close()

    #output weights to a file
    calibrationData={}
    calibrationData['IntegrationRanges'] = [ {'first':fLayer, 'last':lLayer} for fLayer,lLayer in integRanges ]
    calibrationData['IntegrationSubDet'] = [ item for item in integSD ]
    calibrationData['MIPinKeV']          = [ item for item in mipEn ]
    calibrationData['DefaultWeights']    = [ item for item in defWeights ]
    calibrationData['OptimWeights']      = [ item for item in optimWeights ]
    with io.open('%s/weights.dat'%outDir, 'w', encoding='utf-8') as f: f.write(unicode(json.dumps(calibrationData, sort_keys = True, ensure_ascii=False, indent=4)))

    #prepare energy estimators
    funcArgs='{edep0'
    rawEnFunc,weightEnFunc,optimWeightEnFunc='edep0','%f*edep0'%defWeights[0],'%f*edep0'%optimWeights[0]
    for ireg in xrange(1,len(integRanges)) :
        funcArgs          += ',edep%d'%(ireg)
        rawEnFunc         += '+edep%d'%(ireg)
        weightEnFunc      += '+%f*edep%d'%(defWeights[ireg],ireg)
        optimWeightEnFunc += '+%f*edep%d'%(optimWeights[ireg],ireg)
    weightEnFunc      += '+%f'%defWeights[len(defWeights)-1]
    optimWeightEnFunc += '+%f'%optimWeights[len(optimWeights)-1]
    funcArgs=funcArgs+',en}'
    vars=[
        [ws.data('data').addColumn(ws.factory("RooFormulaVar::rawEnFunc('("+rawEnFunc+")/en',"+funcArgs+")" )),                '#Sigma E_{i}'],
        [ws.data('data').addColumn(ws.factory("RooFormulaVar::weightEnFunc('("+weightEnFunc+")/en',"+funcArgs+")")),           '#Sigma w_{i}E_{i}' ],
        [ws.data('data').addColumn(ws.factory("RooFormulaVar::optimWeightEnFunc('("+optimWeightEnFunc+")/en',"+funcArgs+")")), '#Sigma w^{optim}_{i}E_{i}']
        ]
    enEstimatorsSet=ROOT.RooArgSet(ws.var('pt'),ws.var('eta'),ws.var('en'))
    for v in vars:
        vName=v[0].GetName().replace('Func','')
        enEstimatorsSet.add( ws.factory('%s[0,0,999999999]'%vName) )
    getattr(ws,'import')( ROOT.RooDataSet('fitdata','fitdata',enEstimatorsSet) )

    #create the fit dataset (skip direct usage of RooFormulaVar in a fit) and store final value
    for ientry in xrange(0,ws.data('data').numEntries()):
        entryVars=ws.data('data').get(ientry)
        for baseVar in ['pt','eta','en']: ws.var(baseVar).setVal( entryVars.find(baseVar).getVal() )
        newEntry=ROOT.RooArgSet( ws.var('pt'),  ws.var('eta'), ws.var('en') )
        for v in vars:
            vName=v[0].GetName().replace('Func','')
            ws.var(vName).setVal( v[0].getVal() )
            newEntry.add( ws.var(vName) )
        ws.data('fitdata').add( newEntry )
    
    #prepare to fit the energy slices
    calibFunc=[]
    resFunc=[]
    c=ROOT.TCanvas('c','c',1400,800)
    varCtr=0
    for var,varTitle in vars:
        varCtr+=1
        vName=var.GetName().replace('Func','')

        #run resolution fits in different rapidity ranges
        for iyrang in xrange(1,len(pGunYRanges)-2):

            yRang=pGunYRanges[iyrang]

            nv=len(calibFunc)
            calibFunc.append(ROOT.TGraphErrors())
            calibFunc[nv].SetMarkerStyle(20+iyrang) #nv+4*nv%2)
            calibFunc[nv].SetTitle('%s |#eta|=%3.1f'%(varTitle,yRang[0]))
            calibFunc[nv].SetFillStyle(0)
            calibFunc[nv].SetMarkerColor(varCtr+1)
            calibFunc[nv].SetLineColor(varCtr+1)
            calibFunc[nv].SetName('calib_%s_%d'%(vName,iyrang))
            resFunc.append(calibFunc[nv].Clone('resol_%s_%d'%(vName,iyrang)))

            c.Clear()
            c.Divide(len(pGunPtRanges)/3+1,2)
            ipad=0
            for ptRang in pGunPtRanges:
                ipad=ipad+1
                postfix='fit%d%d%d'%(ipad,nv,iyrang)

                #cut data
                redData=ws.data('fitdata').reduce('pt>=%f && pt<=%f && eta>=%f && eta<=%f'%(
                    ptRang[0]-ptRang[1], ptRang[0]+ptRang[1],
                    yRang[0]-yRang[1],   yRang[0]+yRang[1]))
                if redData.numEntries()<10 :
                    ipad=ipad-1
                    continue
                    
                #generator level information
                fitDataMean,   fitDataSigma   = redData.mean(ws.var(vName)), redData.sigma(ws.var(vName))
                fitDataEnMean, fitDataEnSigma = redData.mean(ws.var('en')),     redData.sigma(ws.var('en'))

                #define PDF and ranges to fit/show
                ws.var(vName).setRange(postfix,fitDataMean-4*fitDataSigma,fitDataMean+4*fitDataSigma)
                ws.var(vName).setRange('fit'+postfix,fitDataMean-2*fitDataSigma,fitDataMean+2*fitDataSigma)
                ws.factory('Gaussian::g_%s(%s,mean_%s[%f,%f,%f],sigma_%s[%f,%f,%f])'%
                           (postfix,vName,
                            postfix,fitDataMean,fitDataMean-5*fitDataSigma,fitDataMean+5*fitDataSigma,
                            postfix,fitDataSigma,fitDataSigma/2,fitDataSigma*2)
                    )

                #fit
                fres=ws.pdf('g_%s'%postfix).fitTo( redData, ROOT.RooFit.Range('fit'+postfix), ROOT.RooFit.Save(True) )
                meanFit       = ws.var('mean_%s'%postfix).getVal()
                meanFitError  = ws.var('mean_%s'%postfix).getError()
                sigmaFit      = ws.var('sigma_%s'%postfix).getVal()
                sigmaFitError = ws.var('sigma_%s'%postfix).getError()
                if meanFit<0:
                    ipad=ipad-1
                    continue

                #save results
                np=calibFunc[nv].GetN()
                calibFunc[nv].SetPoint(np,fitDataEnMean,meanFit*fitDataEnMean)
                calibFunc[nv].SetPointError(np,fitDataEnSigma,meanFitError*fitDataEnMean)
                resFunc[nv].SetPoint(np,1/ROOT.TMath.Sqrt(fitDataEnMean),sigmaFit/meanFit)
                resFunc[nv].SetPointError(np,fitDataEnSigma/(2*ROOT.TMath.Power(fitDataEnMean,1.5)),sigmaFitError/meanFit)
            
                #show the result
                p=c.cd(ipad)
                frame=ws.var(vName).frame(ROOT.RooFit.Range(postfix))
                redData.plotOn(frame)
                ws.pdf('g_%s'%postfix).plotOn(frame,ROOT.RooFit.Range(postfix))
                frame.Draw()
                frame.GetXaxis().SetNdivisions(5)
                frame.GetXaxis().SetTitle(varTitle + ' / Energy')
                frame.GetYaxis().SetTitle('Events')
                frame.GetYaxis().SetTitleOffset(1.4)
                frame.GetYaxis().SetRangeUser(0,2*frame.GetMaximum())
                pt=MyPaveText('[E_{beam}=%3.1f GeV, |#eta|=%3.1f]\\<E>=%3.1f RMS=%3.2f\\#mu=%3.1f #sigma=%3.2f'%(fitDataEnMean,yRang[0],fitDataMean,fitDataSigma,meanFit,sigmaFit),
                              0.18,0.9,0.5,0.6)
                pt.SetTextFont(42)
                pt.SetTextSize(0.06)
                if ipad==1:
                    pt=MyPaveText('CMS simulation')
                    pt.SetTextSize(0.08)

            for ext in ['png','pdf','C'] : c.SaveAs('%s/efits_%s_eta%3.1f.%s'%(outDir,vName,yRang[0],ext))
            
            
    return showEMcalibrationResults(calibFunc=calibFunc,resFunc=resFunc,outDir=outDir)

"""
steer 
"""
def main(argv=None):

    customROOTstyle()
    #ROOT.gROOT.SetBatch(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    #genM=0.000511
    #url='cmssw/SingleElectron_SLHC13_30um_SimHits.root'

    genM=0.1396
    url='cmssw/SinglePion_SLHC13_30um_SimHits.root'

    runResolutionStudy(url=url,genM=genM)
    
if __name__ == "__main__":
    sys.exit(main())
