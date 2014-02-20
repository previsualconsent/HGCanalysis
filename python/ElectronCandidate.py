from UserCode.HGCanalysis.PlotUtils import *

"""
Wrapper for electron energy deposits, etc.
"""
class ElectronCandidate:

    def __init__(self):
        self.edeps=[]
        self.edeps_etaphi={}
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

    def setGeneratorLevelInfo(self,gen_en,gen_pt,gen_eta):
        self.gen_en=gen_en
        self.gen_pt=gen_pt
        self.gen_eta=gen_eta

    def setLongitudinalProfile(self,edeps,localOverburden):
        self.edeps=edeps
        self.localOverburden=localOverburden
        
    def setTransverseProfile(self,edeps_etaphi):
        self.edeps_etaphi=edeps_etaphi
        
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
            
            pt=MyPaveText('E^{gen}=%3.0f GeV\\p_{T}^{gen}=%3.0f GeV, #eta^{gen}=%3.0f\\#Sigma E_{i}=%3.0f\\#Sigma w_{i}E_{i}=%3.0f\\Fit E=%3.0f'
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

        #save some plots
        c=ROOT.TCanvas("cetaphi","cetaphi",1000,1000);
        c.Divide(3,3)
        sumEtaPhi=ROOT.TH2F('sumetaphi',';Pseudo-rapidity; #phi [rad]; Energy',50,-0.5,0.5,300,-1.5,1.5)
        sumEtaPhi.GetXaxis().SetLabelSize(0.08);
        sumEtaPhi.GetXaxis().SetTitleSize(0.08);
        sumEtaPhi.GetXaxis().SetTitleOffset(1.0);
        sumEtaPhi.GetXaxis().SetNdivisions(5);
        sumEtaPhi.GetYaxis().SetLabelSize(0.08);
        sumEtaPhi.GetYaxis().SetTitleSize(0.08);
        sumEtaPhi.GetYaxis().SetTitleOffset(1.2);
        sumEtaPhi.GetYaxis().SetNdivisions(5);
        sumEtaPhi.GetZaxis().SetLabelSize(0.07);
        sumEtaPhi.GetZaxis().SetTitleSize(0.08);
        sumEtaPhi.GetZaxis().SetTitleOffset(1.0);
        sumEtaPhi.GetZaxis().SetNdivisions(5);
        sumEtaPhi.GetZaxis().SetRangeUser(1e2,1e8);
        sumEtaPhi_perLayer=[]
        layersToSample=[3,6,8,12,14,16,20,24]
        layersToSample3D=[8,14,20]
        ipad=0
        for layer in self.edeps_etaphi:

            nhistos=len(sumEtaPhi_perLayer)
            sumEtaPhi_perLayer.append(sumEtaPhi.Clone('sumetaphi_'+str(layer)))
            sumEtaPhi_perLayer[nhistos].Reset('ICE')
            for edep in self.edeps_etaphi[layer]:
                print edep[0],edep[1],edep[2]
                sumEtaPhi_perLayer[nhistos].Fill(edep[0],edep[1],edep[2])
            weight=self.localOverburden[layer-1]/self.localOverburden[0]
            sumEtaPhi.Add( sumEtaPhi_perLayer[ nhistos ], weight )

            if layer in layersToSample:
                ipad=ipad+1
                p=c.cd(ipad)
                p.SetLogz();
                p.SetBottomMargin(0.2);
                p.SetRightMargin(0.05);
                p.SetLeftMargin(0.2);
                p.SetTopMargin(0.05);
                if layer in layersToSample3D:
                    ROOT.gStyle.SetPalette(55);
                    #sumEtaPhi_perLayer[ nhistos ].Draw("surf2fb 0");
                    sumEtaPhi_perLayer[ nhistos ].Draw("lego2fb 0");
                    sumEtaPhi_perLayer[ nhistos ].GetZaxis().SetRangeUser(1,5e3);
                    pt=MyPaveText('[Layer %d]'%layer,0.2,0.75,0.6,0.8)
                    pt.SetTextFont(42)
                    pt.SetTextAlign(12)
                    pt.SetTextColor(16)
                    pt.SetTextSize(0.09)
                else :
                    opt='col'
                    if ipad==1 : opt+='z'
                    sumEtaPhi_perLayer[ nhistos ].Draw(opt)
                    sumEtaPhi_perLayer[ nhistos ].GetZaxis().SetRangeUser(1,5e3);
                    pt=MyPaveText('[Layer %d]'%layer,0.2,0.25,0.6,0.3)
                    pt.SetTextFont(42)
                    pt.SetTextAlign(12)
                    pt.SetTextColor(16)
                    pt.SetTextSize(0.09)
                if ipad==1:
                    p.SetRightMargin(0.2);
                    sumEtaPhi_perLayer[ nhistos ].GetZaxis().SetLabelSize(0.07);
                    sumEtaPhi_perLayer[ nhistos ].GetZaxis().SetTitleSize(0.08);
                    sumEtaPhi_perLayer[ nhistos ].GetZaxis().SetTitleOffset(0.9);
                if ipad==1:
                    MyPaveText('CMS simulation\\E^{gen}=%3.0f GeV\\ p_{T}^{gen}=%3.0f GeV, #eta^{gen}=%3.0f'
                               %(self.gen_en,self.gen_pt,self.gen_eta),0.22,0.8)

        # show the sum in the last pad
        p=c.cd(9)
        p.SetLogz();
        p.SetBottomMargin(0.2);
        p.SetRightMargin(0.05);
        p.SetLeftMargin(0.2);
        p.SetTopMargin(0.05);
        sumEtaPhi.Draw("lego2fb 0");
        sumEtaPhi.GetZaxis().SetRangeUser(1,1e4)
        pt=MyPaveText('[#Sigma w_{i} E_{i}]')
        pt.SetTextFont(42)
        pt.SetTextAlign(12)
        pt.SetTextColor(16)
        pt.SetTextSize(0.09)

        c.Modified()
        c.Update()
        title='cmssw_pt%d_eta%d_etaphi'%(self.gen_pt,self.gen_eta)
        raw_input()
        c.SaveAs(title+'.png')
        sumEtaPhi.Delete()
        for h in sumEtaPhi_perLayer: h.Delete()
                    
