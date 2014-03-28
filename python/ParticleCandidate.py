from UserCode.HGCanalysis.PlotUtils import *

"""
Wrapper for electron energy deposits, etc.
"""
class ParticleCandidate:

    def __init__(self):
        self.edeps=[]
        self.edeps_xy={}
        self.localOverburden=[]
        self.gen_en=-1
        self.gen_pt=-1
        self.gen_eta=-1
        self.totalEn=0
        self.totalEnVsRho=[0,0,0,0]
        self.totalRawEn=0
        self.totalFitEn=0
        self.showerFunc=None
        self.chi2=0
        self.ndof=0
        self.showerMax=0

    def setGeneratorLevelInfo(self,gen_en,gen_pt,gen_eta,gen_phi):
        self.gen_en=gen_en
        self.gen_pt=gen_pt
        self.gen_eta=gen_eta
        self.gen_phi=gen_phi

    def setLongitudinalProfile(self,edeps,localOverburden):
        self.edeps=edeps
        self.localOverburden=localOverburden
        
    def setTransverseProfile(self,edeps_xy):
        self.edeps_xy=edeps_xy
        
    def getEnergyAt(self,layer,applyCorr=False):
        if not applyCorr : return self.edeps[layer]
        return self.edeps[layer]*self.localOverburden[layer]/self.localOverburden[0]

    def buildLongitudinalProfile(self,drawCanvas=False):

        #build the longitudinal graph profile and fit it
        self.longProfile=ROOT.TGraphErrors()
        self.longProfile.SetName('longprof')
        self.totalEn=0
        self.totalRawEn=0
        totalX0=0
        for layer in xrange(0,len(self.edeps)) :
            localX0          = self.localOverburden[layer]
            totalX0         += localX0
            rawEn            = self.getEnergyAt(layer,False)
            self.totalEn    += self.getEnergyAt(layer,True)
            self.totalRawEn += rawEn

            np=self.longProfile.GetN()
            self.longProfile.SetPoint(np,totalX0,rawEn)
            self.longProfile.SetPointError(np,0,ROOT.TMath.Sqrt(rawEn))

        #fit a Gamma-like function to the shower profile
        self.showerFunc=ROOT.TF1('showerfunc','[0]*pow(x,[1])*exp(-[2]*x)',0,totalX0)
        self.showerFunc.SetParLimits(1,0.,100.);
        self.showerFunc.SetParLimits(2,0.,100.);
        self.showerFunc.SetLineColor(3);
        self.longProfile.Fit(self.showerFunc,"RQ+"); 
        self.chi2=self.showerFunc.GetChisquare()
        self.ndof=self.showerFunc.GetNDF()
        if self.showerFunc.GetParameter(2) > 0: 
            self.showerMax=self.showerFunc.GetParameter(1)/self.showerFunc.GetParameter(2)
        else:
            self.showerMax=0
        self.totalFitEn=self.showerFunc.Integral(0,totalX0)

        #show in canvas if required
        if drawCanvas :
            
            c=ROOT.TCanvas('c','c',500,500)
            c.SetTopMargin(0.05)
            self.longProfile.Draw('ap')
            self.longProfile.SetMarkerStyle(20)
            self.longProfile.GetXaxis().SetTitle('Transversed thickness [1/X_{0}^{eff}]')
            self.longProfile.GetYaxis().SetTitle('Energy')
            self.longProfile.GetXaxis().SetLabelSize(0.04)
            self.longProfile.GetYaxis().SetLabelSize(0.04)
            self.longProfile.GetXaxis().SetTitleSize(0.05)
            self.longProfile.GetYaxis().SetTitleSize(0.05)
            self.longProfile.GetYaxis().SetTitleOffset(1.4)
            
            MyPaveText('CMS simulation')
            
            pt=MyPaveText('E^{gen}=%3.0f GeV\\p_{T}^{gen}=%3.0f GeV, #eta^{gen}=%3.1f\\#Sigma E_{i}=%3.0f\\#Sigma w_{i}E_{i}=%3.0f\\Fit E=%3.0f'
                          %(self.gen_en,self.gen_pt,self.gen_eta,self.totalRawEn,self.totalEn,self.totalFitEn),
                          0.6,0.6,0.9,0.9)
            pt.AddText('Shower max=%3.1f'%self.showerMax)
            pt.AddText('#chi^{2}/ndof=%3.0f/%d'%(self.chi2,self.ndof)) 
            pt.SetTextFont(42)
            pt.SetTextSize(0.03)

            c.Modified()
            c.Update()
            title='cmssw_pt%d_eta%d'%(self.gen_pt,self.gen_eta)
            c.SaveAs(title+'.png')
            c.SaveAs(title+'.pdf')
        
        #all done here
        return self.longProfile

    def buildTransverseProfile(self,debug):

        #for layer in self.edeps_etaphi:
        #    weight=self.localOverburden[layer-1]/self.localOverburden[0]
        #    for edep in self.edeps_etaphi[layer]:
        #        rho=ROOT.TMath.Sqrt(ROOT.TMath.Power(edep[1],2)+ROOT.TMath.Power(edep[2],2))
        #        val=edep[0]*weight

        #self.totalEnVsRho[0] += val
        #        if rho<10: self.totalEnVsRho[1] += val
        #        if rho<20: self.totalEnVsRho[2] += val
        #        if rho<1.70*ROOT.TMath.Exp(0.106*layer) : self.totalEnVsRho[3] += val

                
        if not debug : return



        #can probably be done with numpy...
        avgX,avgY,totalEn=0,0,0
        for layer in self.edeps_xy:
            for edep in self.edeps_xy[layer]:
                if not edep[5]: continue
                totalEn=totalEn+edep[2]
                dR=ROOT.TMath.Sqrt(edep[0]*edep[0]+edep[1]*edep[1])
                avgX=avgX+edep[3]*edep[2]
                avgY=avgY+edep[4]*edep[2]
        if totalEn>0:
            avgX=avgX/totalEn
            avgY=avgY/totalEn



        #save some plots
        c=ROOT.TCanvas("cxy","cxy",1000,1000);
        c.Divide(3,3)
        sumXY=ROOT.TH2F('sumxy',';x [mm]; y [mm]; Energy/MIP',40,avgX-295,avgX+195,40,avgY-205,avgY+195)
        sumXY.GetXaxis().SetLabelSize(0.06);
        sumXY.GetXaxis().SetTitleSize(0.07);
        sumXY.GetXaxis().SetTitleOffset(1.0);
        sumXY.GetXaxis().SetNdivisions(5);
        sumXY.GetYaxis().SetLabelSize(0.06);
        sumXY.GetYaxis().SetTitleSize(0.07);
        sumXY.GetYaxis().SetTitleOffset(1.2);
        sumXY.GetYaxis().SetNdivisions(5);
        sumXY.GetZaxis().SetLabelSize(0.06);
        sumXY.GetZaxis().SetTitleSize(0.07);
        sumXY.GetZaxis().SetTitleOffset(1.0);
        sumXY.GetZaxis().SetNdivisions(5);
        sumXY.GetZaxis().SetRangeUser(1e2,1e8);
        sumXY_perLayer={}
        for layer in self.edeps_xy:
            print layer,len(self.edeps_xy[layer])
            sumXY_perLayer[layer]=sumXY.Clone('sumxy_'+str(layer))
            sumXY_perLayer[layer].Reset('ICE')
            for edep in self.edeps_xy[layer]:
                sumXY_perLayer[layer].Fill(edep[3],edep[4],edep[2])
            weight=1.0 #self.localOverburden[layer-1]/self.localOverburden[0]
            sumXY.Add( sumXY_perLayer[ layer ], weight )

        layersToSample=[3,6,8,12,14,16,20,24]
        layersToSample3D=[8,14,20]
        ipad=0
        sumXY_toDisplay={}
        maxEn=0
        print len(sumXY_perLayer)
        for layer in layersToSample:
            sumXY_toDisplay[layer]=sumXY_perLayer[layer*3].Clone('sumxyint_%d'%layer)
            sumXY_toDisplay[layer].Add(sumXY_perLayer[layer*3+1])
            #sumXY_toDisplay[layer].Add(sumXY_perLayer[layer*3+2])
            maxEn=max(maxEn,sumXY_toDisplay[layer].GetMaximum())


        for layer in sumXY_toDisplay:
            print layer
            if layer in layersToSample:
                ipad=ipad+1
                p=c.cd(ipad)
                #p.SetLogz();
                p.SetBottomMargin(0.2);
                p.SetRightMargin(0.05);
                p.SetLeftMargin(0.2);
                p.SetTopMargin(0.05);
                if layer in layersToSample3D:
                    ROOT.gStyle.SetPalette(55);
                    #sumXY_perLayer[ nhistos ].Draw("surf2fb 0");
                    sumXY_toDisplay[ layer ].Draw("lego2fb 0");
                    sumXY_toDisplay[ layer ].GetZaxis().SetRangeUser(0,maxEn);
                    pt=MyPaveText('[Layer %d]'%layer,0.2,0.75,0.6,0.8)
                    pt.SetTextFont(42)
                    pt.SetTextAlign(12)
                    #pt.SetTextColor(16)
                    pt.SetTextSize(0.07)
                else :
                    opt='col'
                    if ipad==1 : opt+='z'
                    sumXY_toDisplay[ layer ].Draw(opt)
                    sumXY_toDisplay[ layer ].GetZaxis().SetRangeUser(0,maxEn);
                    pt=MyPaveText('[Layer %d]'%layer,0.2,0.25,0.6,0.3)
                    pt.SetTextFont(42)
                    pt.SetTextAlign(12)
                    #pt.SetTextColor(16)
                    pt.SetTextSize(0.07)
                if ipad==1:
                    p.SetRightMargin(0.2);
                    sumXY_toDisplay[ layer ].GetZaxis().SetLabelSize(0.07);
                    sumXY_toDisplay[ layer ].GetZaxis().SetTitleSize(0.08);
                    sumXY_toDisplay[ layer ].GetZaxis().SetTitleOffset(0.9);
                if ipad==1:
                    MyPaveText('CMS simulation, #sqrt{s}=14 TeV',0.12).SetTextSize(0.05)
                    
        # show the sum in the last pad
        p=c.cd(9)
        #p.SetLogz();
        p.SetBottomMargin(0.2);
        p.SetRightMargin(0.05);
        p.SetLeftMargin(0.2);
        p.SetTopMargin(0.05);
        sumXY.Draw("lego2fb 0");
        sumXY.GetZaxis().SetRangeUser(0,maxEn*2)
        pt=MyPaveText('[#Sigma E_{i}]',0.15,0.8)
        pt.SetTextFont(42)
        pt.SetTextAlign(12)
        pt.SetTextColor(16)
        pt.SetTextSize(0.08)
        pt=MyPaveText('E^{gen}=%3.0f GeV, p_{T}^{gen}=%3.0f GeV\\ #eta^{gen}=%3.1f, #phi^{gen}=%3.1f'
                      %(self.gen_en,self.gen_pt,self.gen_eta,self.gen_phi),0.4,0.8)
        pt.SetTextSize(0.05)
        pt.SetTextFont(42)

        c.Modified()
        c.Update()
        title='cmssw_pt%d_eta%d_xy'%(self.gen_pt,self.gen_eta)
        raw_input()
        c.SaveAs(title+'.png')
        sumXY.Delete()
        #for h in sumXY_perLayer: h.Delete()
                    
