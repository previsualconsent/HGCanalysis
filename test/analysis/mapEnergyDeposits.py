#!/usr/bin/env python

import ROOT
import sys
from ROOT import *
from UserCode.HGCanalysis.PlotUtils import *

customROOTstyle()
gROOT.SetBatch(False)#True)
gStyle.SetPalette(55)

inF=TFile.Open('HGCSimHitsAnalysis.root')
Events=inF.Get('hgcSimHitsAnalyzer/HGC')

#generate the event sets
totalEvts=Events.GetEntriesFast()
nPu=int(sys.argv[2])
nSets=totalEvts
curEv=1
evRanges=[]
while True:
    iPu=gRandom.Poisson(nPu)
    evRanges.append([curEv,curEv+iPu])
    curEv=curEv+iPu+1
    if curEv>=totalEvts or curEv>2500:
        break
print 'Generating %d event sets'%len(evRanges)

#quantiles to profile at each layer for the pileup event sets
import array
xq=array.array('d',[0.68,0.95])
yq=array.array('d',[0,00,0.00])

#get number of sectors per layer
sectorPerLayer={}
for i in xrange(0,1):
    Events.GetEntry(i)
    for n in xrange(0,Events.nee):
        try:
            sectorPerLayer[Events.ee_layer[n]].append(Events.ee_sec[n])
        except:
            sectorPerLayer[Events.ee_layer[n]]=[Events.ee_sec[n]]
            
#yranges=[[-500,-495],[-250,-245],[0,5],[245,250],[495,500]]
yranges=[[-500,-495],[-250,-245],[0,5],[495,500]]

#iterate over layers
for l in sectorPerLayer:

    if l!=int(sys.argv[1]) : continue
            
    sectorPerLayer[l]=list(set(sectorPerLayer[l]))
    nSectors=len(sectorPerLayer[l])
    print '\t Analyzing layer %d over %d sectors'%(l,nSectors)

    c=TCanvas('c','c',500,500)

    var='ee_y:ee_subsec*ee_x >> hlayer(500,-500,500,1000,-1000,1000)'
    cut='(ee_edep*1e6/56)*(1./(%f*%f))*(ee_layer==%d && event<%d)'%(len(evRanges),nSectors,l,evRanges[len(evRanges)-1][1])
    Events.Draw(var,cut,'colz')
    averageHitEn=c.GetPrimitive('hlayer').Clone('averageHitEn')
    averageHitEn.SetDirectory(0)

    hitDist=TH2F('hitdist',';#hits;y[mm];Events',100,0,100,len(yranges),0,len(yranges))
    hitEnDist=TH2F('hitendist',';Energy / 56 keV;y[mm];<Hits>',200,0,100,len(yranges),0,len(yranges))
    q50Dist=TH2F('q50en',';median energy / 56 keV;y[mm];Events',100,0,100,len(yranges),0,len(yranges))
    q95Dist=TH2F('q95en',';95% quantile energy / 56 keV;y[mm];Events',100,0,100,len(yranges),0,len(yranges))    
    for iy in xrange(0,len(yranges)):

        print ' \t Starting with ',yranges[iy]

        #if iy>0: continue
        
        y=yranges[iy]
        hitDist.GetYaxis().SetBinLabel(iy+1,'[%d,%d]'%(y[0],y[1]))
        hitEnDist.GetYaxis().SetBinLabel(iy+1,'[%d,%d]'%(y[0],y[1]))
        q50Dist.GetYaxis().SetBinLabel(iy+1,'[%d,%d]'%(y[0],y[1]))
        q95Dist.GetYaxis().SetBinLabel(iy+1,'[%d,%d]'%(y[0],y[1]))

        #overlay minbias events for each y strip
        for iSet in evRanges:

            cut='ee_edep*1e6>56 && ee_layer==%d && event>=%d && event<%d && ee_y>=%f && ee_y<%f'%(l,iSet[0],iSet[1],y[0],y[1])
            nhits=Events.GetEntries(cut)
            hitDist.Fill(nhits,iy)
            
            var='ee_edep*1e6/56 >> hen(200,0,100)'
            opt='hist'
            Events.Draw(var,cut,opt)        
            h=c.GetPrimitive('hen')

            for xbin in xrange(1,hitEnDist.GetXaxis().GetNbins()) :
                hitEnDist.Fill(h.GetBinCenter(xbin),iy, h.GetBinContent(xbin)/len(evRanges))

            #compute the quantiles for this set of pu
            h.GetQuantiles(len(xq),yq,xq)
            q50Dist.Fill(yq[0],iy)
            q95Dist.Fill(yq[1],iy)


    #show results for this layer
    cl=TCanvas('cl','cl',1000,1000)
    cl.Divide(2,2)

    p=cl.cd(1)
    p.SetRightMargin(0.15)
    p.SetLogz()
    averageHitEn.Draw('colz')
    averageHitEn.GetXaxis().SetNdivisions(5)
    averageHitEn.GetYaxis().SetNdivisions(5)
    averageHitEn.GetZaxis().SetNdivisions(5)
    averageHitEn.GetZaxis().SetTitle('Energy/(56 keV)')
    averageHitEn.GetYaxis().SetTitle('y [mm]')
    averageHitEn.GetXaxis().SetTitle('x [mm]')
    averageHitEn.GetZaxis().SetTitleSize(0.05)
    averageHitEn.GetYaxis().SetTitleSize(0.05)
    averageHitEn.GetXaxis().SetTitleSize(0.05)
    averageHitEn.GetZaxis().SetLabelSize(0.04)
    averageHitEn.GetYaxis().SetLabelSize(0.04)
    averageHitEn.GetXaxis().SetLabelSize(0.04)

    MyPaveText('#splitline{CMS simulation, #sqrt{s}=13 TeV}{Layer %d, <PU>=%d}'%(l,nPu),  0.15,0.92,0.6,0.96)     

    p=cl.cd(2)
    p.SetRightMargin(0.15)
    hitDist.Draw('colz')
    hitDist.GetZaxis().SetTitleSize(0.05)
    hitDist.GetYaxis().SetTitleSize(0.05)
    hitDist.GetXaxis().SetTitleSize(0.05)
    hitDist.GetZaxis().SetLabelSize(0.04)
    hitDist.GetYaxis().SetLabelSize(0.04)
    hitDist.GetXaxis().SetLabelSize(0.04)

    p=cl.cd(3)
    p.SetRightMargin(0.15)
    hitEnDist.Draw('colz')
    hitEnDist.GetXaxis().SetRangeUser(0,20)
    hitEnDist.GetZaxis().SetTitleSize(0.05)
    hitEnDist.GetYaxis().SetTitleSize(0.05)
    hitEnDist.GetXaxis().SetTitleSize(0.05)
    hitEnDist.GetZaxis().SetLabelSize(0.04)
    hitEnDist.GetYaxis().SetLabelSize(0.04)
    hitEnDist.GetXaxis().SetLabelSize(0.04)

    p=cl.cd(4)
    p.Divide(1,2)
    sp=p.cd(1)
    sp.SetRightMargin(0.15)
    q50Dist.Draw('colz')
    q50Dist.GetZaxis().SetTitleSize(0.05)
    q50Dist.GetYaxis().SetTitleSize(0.05)
    q50Dist.GetXaxis().SetTitleSize(0.05)
    q50Dist.GetZaxis().SetLabelSize(0.04)
    q50Dist.GetYaxis().SetLabelSize(0.04)
    q50Dist.GetXaxis().SetLabelSize(0.04)
    
    sp=p.cd(2)
    sp.SetRightMargin(0.15)

    q95Dist.Draw('colz')
    q95Dist.GetZaxis().SetTitleSize(0.06)
    q95Dist.GetYaxis().SetTitleSize(0.06)
    q95Dist.GetXaxis().SetTitleSize(0.06)
    q95Dist.GetZaxis().SetLabelSize(0.05)
    q95Dist.GetYaxis().SetLabelSize(0.05)
    q95Dist.GetXaxis().SetLabelSize(0.05)

    cl.cd()
    cl.Modified()
    cl.Update()
    cl.SaveAs('layer_%d_pu_%d.png'%(l,nPu))
    cl.SaveAs('layer_%d_pu_%d.root'%(l,nPu))

    raw_input()

inF.Close()


