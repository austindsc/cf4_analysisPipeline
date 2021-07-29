#!/usr/local/bin/python3

#stock python stuff
import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
import time as chrono

## this is pyroot
from ROOT import ROOT
from ROOT import TCanvas, TGraph, gApplication, TPad, TRatioPlot, TLine, TText, TLatex, gStyle, TLegend
from ROOT import TGraphErrors, TF1, TFile
from ROOT import gROOT
from ROOT import TH1

import sys
# sys.path.append('..')
from SBCcode.DataHandling import ReadBinary as rb

#### compare SPE charges for with/without source ###

P2DATA = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/atHome/processed/"

f_s55bkg = TFile.Open(P2DATA + "1123_0.root")

b1 = (f_s55bkg.Get("hq_flat")).Clone()
b2 = (f_s55bkg.Get("hq_Earlyflat")).Clone()

b3 = b1.Clone()
b3.Add(b2)

f_s55Cs = TFile.Open(P2DATA + "1123_1.root")
c1 = (f_s55Cs.Get("hq_flat")).Clone()
c2 = (f_s55Cs.Get("hq_Earlyflat")).Clone()

c3 = b1.Clone()
c3.Add(c2)

c1 = TCanvas( 'c1', 'c1', 600, 600)

b1.SetLineColor(2)
b1.GetXaxis().SetRangeUser(1, b3.GetMaximum())
b1.Draw("hist")

b2.SetLineColor(3)
b2.Draw("hist same")

b3.Draw("hist same")

c3.SetLineColor(6)
c3.Draw("hist same")

c1.Modified()
c1.Update()
