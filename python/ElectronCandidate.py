from UserCode.HGCanalysis.PlotUtils import *

"""
Wrapper for electron energy deposits, etc.
"""
class ElectronCandidate:
    def __init__(self):
        self.edeps=[]
        self.localOverburden=[]
        self.gen_en=-1
        self.gen_pt=-1
        self.gen_eta=-1
        self.totalEn=0
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
    def getEnergyAt(self,layer,applyCorr=False):
        if not applyCorr : return self.edeps[layer]
        return self.edeps[layer]*self.localOverburden[layer]
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
        self.showerMax=self.showerFunc.GetParameter(1)/self.showerFunc.GetParameter(2)
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

            title='cmssw_pt%d_eta%d'%(self.gen_pt,self.gen_eta)
            c.SaveAs(title+'.png')
            c.SaveAs(title+'.pdf')
        
        #all done here
        return self.longProfile
