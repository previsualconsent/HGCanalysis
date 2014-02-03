import ROOT

"""

"""
class MyPaveText(ROOT.TPaveText):
    def __init__(self,title,x1=0.1,y1=0.95,x2=0.6,y2=0.99):
        ROOT.TPaveText.__init__(self,x1,y1,x2,y2,'brNDC')
        self.SetFillStyle(0)
        self.SetBorderSize(0)
        self.SetTextAlign(12)
        self.SetTextSize(0.04)
        for t in title.split('\\') : self.AddText(t)
        self.Draw()


"""
Style options mostly from CMS's tdrStyle.C
"""
def customROOTstyle() :
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(False)
    ROOT.gStyle.SetPadTopMargin(0.1);
    ROOT.gStyle.SetPadBottomMargin(0.13);
    ROOT.gStyle.SetPadLeftMargin(0.15);
    ROOT.gStyle.SetPadRightMargin(0.02);
    ROOT.gStyle.SetLabelColor(1, "XYZ");
    ROOT.gStyle.SetLabelFont(42, "XYZ");
    ROOT.gStyle.SetLabelOffset(0.007, "XYZ");
    ROOT.gStyle.SetLabelSize(0.05, "XYZ");
    ROOT.gStyle.SetTitleSize(0.05, "XYZ");
    ROOT.gStyle.SetAxisColor(1, "XYZ");
    ROOT.gStyle.SetStripDecimals(True);
    ROOT.gStyle.SetTickLength(0.03, "XYZ");
    ROOT.gStyle.SetNdivisions(510, "XYZ");
    ROOT.gStyle.SetPadTickX(0);
    ROOT.gStyle.SetPadTickY(0);
    ROOT.gStyle.SetMarkerStyle(20);
    ROOT.gStyle.SetHistLineColor(1);
    ROOT.gStyle.SetHistLineStyle(0);
    ROOT.gStyle.SetHistLineWidth(1);
    ROOT.gStyle.SetFrameBorderMode(0);
    ROOT.gStyle.SetFrameBorderSize(1);
    ROOT.gStyle.SetFrameFillColor(0);
    ROOT.gStyle.SetFrameFillStyle(0);
    ROOT.gStyle.SetFrameLineColor(1);
    ROOT.gStyle.SetFrameLineStyle(1);
    ROOT.gStyle.SetFrameLineWidth(1);
