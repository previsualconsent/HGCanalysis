#!/usr/bin/env python

import ROOT
from ROOT import *
from UserCode.HGCanalysis.PlotUtils import *

customROOTstyle()
gROOT.SetBatch(True)
gStyle.SetPalette(55)

minBiasF=TFile.Open('HGCSimHitsAnalysis_minbias.root')
MinBiasEvents=minBiasF.Get('hgcSimHitsAnalyzer/HGC')

inF=TFile.Open('HGCSimHitsAnalysis.root')
Events=inF.Get('hgcSimHitsAnalyzer/HGC')

import sys
event=int(sys.argv[1])

z=1
nmod=360/20

for l in xrange(0,1):

    c=TCanvas('c%d'%l,'c%d'%l,1200,1200)
    if l>0 : MyPaveText('#splitline{Event %d}{Layer %d}'%(event,l),  0.4,0.45,0.6,0.55)
    else   : MyPaveText('#splitline{Event %d}{#Sigma Layer}'%(event),0.4,0.45,0.6,0.55)
    
    c.Divide(9,2)

    for m in xrange(1,nmod+1):

        theta=(m-1)*20*TMath.Pi()/180
        
        xpos='ee_x*ee_subsec'
        ypos='ee_y'
        var='%s:%s'%(ypos,xpos)

        #rotxpos='%s*TMath::Cos(%f)-%s*TMath::Sin(%f)'%(xpos,theta,ypos,theta)
        #rotypos='%s*TMath::Sin(%f)+%s*TMath::Sin(%f)'%(xpos,theta,ypos,theta)
        #var='%s:%s'%(rotypos,rotxpos)
        
        #cut='ee_edep*(ee_layer==%d && ee_zp==%d && ee_module==%d && event==%d)'%(l,z,m,event)
        if l>0 :
            cut  ='ee_edep*(ee_layer==%d && ee_zp==%d && ee_module==%d)'%(l,z,m)
            evcut='ee_edep*(event==%d && ee_layer==%d && ee_zp==%d && ee_module==%d)'%(event,l,z,m)
        else :
            cut  ='(ee_layer<=10 ? 1 : ee_layer<=20 ? 0.8/0.5 : 1.2/0.5)*ee_edep*(ee_zp==%d && ee_module==%d)'%(z,m)
            evcut='(ee_layer<=10 ? 1 : ee_layer<=20 ? 0.8/0.5 : 1.2/0.5)*ee_edep*(event==%d && ee_zp==%d && ee_module==%d)'%(event,z,m)

        pad=c.cd(m)
    
        scaleFactor=4.0/10.0
        baseX=TMath.Cos(theta)*scaleFactor+0.5
        deltaX=-0.06*TMath.Sign(1,TMath.Cos(theta))
        baseY=TMath.Sin(theta)*scaleFactor+0.5
        deltaY=-0.06*TMath.Sign(1,TMath.Sin(theta))
        pad.SetPad('subp%d%d'%(l,m),'Sector=%d'%m,baseX-deltaX,baseY-deltaY,baseX+deltaX,baseY+deltaY,0,0,0)
        pad.SetLogz()
        pad.SetBottomMargin(0)
        pad.SetTopMargin(0.1)
        pad.SetLeftMargin(0)
        pad.SetRightMargin(0.1)
        hname='h%d'%m

        opt='col'
        if m==3 : opt+='z'

        #opt='surf2 FB BB 0'

        MinBiasEvents.Draw(var+'>>%s(500,-500,500,1000,-1000,1000)'%hname,cut,opt)
        Events.Draw(var+'>>+%s'%hname,evcut,opt)
        h=pad.GetPrimitive(hname)
        h.RebinX(10)
        h.RebinY(10)

        tokeV=1e6
        Eeh=3.6e-3
        h.Scale(tokeV/Eeh)
        h.GetXaxis().SetNdivisions(5)
        h.GetYaxis().SetNdivisions(5)
        h.GetZaxis().SetNdivisions(5)
        h.GetZaxis().SetRangeUser(1e3,1e9)
        h.GetZaxis().SetTitle('Energy/E_{eh}')
        h.GetZaxis().SetTitleSize(0.09)
        h.GetZaxis().SetTitleOffset(1.4)
        h.GetZaxis().SetLabelSize(0.09)
        h.GetYaxis().SetLabelSize(0.09)
        h.GetXaxis().SetLabelSize(0.08)
        
        
        txt=MyPaveText('[Sector %d]'%m,0.25,0.8,0.7,0.85)
        txt.SetTextSize(0.09)

    c.Modified()
    c.Update()
    c.SaveAs('evdisplay_%d_%d.png'%(l,event))

minBiasF.Close()
inF.Close()
