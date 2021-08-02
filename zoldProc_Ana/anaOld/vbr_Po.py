#!/usr/local/bin/python3

# system stuff
import sys
sys.path.append("../")
import time

#root stuff
import ROOT
from ROOT import TFile
from ROOT import TCanvas, gApplication, gStyle, TLegend, TPaveText
from ROOT import TF1, TH1F, TGraph, TGraphErrors, TFitResult
from ROOT import TMath, TFormula

import numpy as np

##### Plan ish #####
# 1. get Po events in gas for different voltages
# 2. fit high energy peaks
# 3. plot high energy peak's means
# 4. fit with linear 
# 5. fit with exponential corrected P_extra

## objects ## 
# TH1F	hF_cuts; - charge histograms
# TH1F	hH_cuts; - height histograms
# TH2F	h2_hq_cuts;
# TH2F	h2_qw_cuts;
# TH2F	h2_hw_cuts;

# opening folders
P2DATA = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/CF4_coding/proc1/"

def vbr_Po():
	fold = "gas/"
	n51 = "gasBig_51_Po_1123_x20.root"
	n53 = "gasBig_53_Po_1123_x20.root"
	n55 = "gasBig_55_Po_1123_x20.root"

	f_51 = TFile.Open(P2DATA + fold + n51)
	hq51 = (f_51.Get("hF_cuts")).Clone()

	f_53 = TFile.Open(P2DATA + fold + n53)
	hq53 = (f_53.Get("hF_cuts")).Clone()

	f_55 = TFile.Open(P2DATA + fold + n55)
	hq55 = (f_55.Get("hF_cuts")).Clone()

	# scaling and stuff
	hq51.Rebin(10)
	hq53.Rebin(10)
	hq55.Rebin(10)

	# plotting
	c1 = TCanvas( 'c1', 'c1', 600, 600)

	hq51.GetYaxis().SetRangeUser(1, hq55.GetMaximum()+100)
	hq51.SetLineColor(8)
	hq51.Draw("hist")

	hq53.SetLineColor(9)
	hq53.Draw("hist same")

	hq55.SetLineColor(46)
	hq55.Draw("hist same")

	c1.SetLogy()

	c1.Modified()
	c1.Update()

	ROOT.gApplication.Run()

if __name__ == '__main__':
	vbr_Po()
