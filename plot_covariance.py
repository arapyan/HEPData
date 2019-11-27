import ROOT as rt
import tdr, CMS_lumi
import array
from numpy import *
from ROOT import TCanvas, TGraph, TLine
from ROOT import kBlack, kBlue, kRed, kWhite, kOrange


#set the tdr style
tdr.setTDRStyle()

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
CMS_lumi.lumi_sqrtS = "35.9 fb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPos = 1
if( iPos==0 ): CMS_lumi.relPosX = 0.12


W = 800;
H = 800;

# references for T, B, L, R
T = 0.08*W
B = 0.10*W 
L = 0.14*W
R = 0.10*W

inputfile = rt.TFile("inputs/smp18006/zv_hist.root","READ")

canvas = rt.TCanvas("c2","c2",50,50,W,H)
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetLeftMargin( L/W )
canvas.SetRightMargin( R/W )
canvas.SetTopMargin( T/H )
canvas.SetBottomMargin( B/H )
canvas.SetTickx(1)
canvas.SetTicky(1)

data = inputfile.Get("hcov")

data.GetXaxis().SetTitle("m_{ZV} (GeV)")
data.GetXaxis().SetLabelSize(0.04)
data.GetXaxis().SetTitleSize(0.045)
data.GetXaxis().SetTitleOffset(0.95)
data.GetYaxis().SetTitleOffset(1.4);
data.GetYaxis().SetTitle("m_{ZV} (GeV)")
data.GetYaxis().SetLabelSize(0.04)
data.GetYaxis().SetTitleSize(0.045)

data.Draw("COLZ")
CMS_lumi.CMS_lumi(canvas, 0, iPos)


canvas.cd()
canvas.Update()
canvas.RedrawAxis()

raw_input("Press Enter to end")
