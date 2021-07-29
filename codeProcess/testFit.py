#!/usr/local/bin/python3

import numpy as np
import os
import matplotlib.pyplot as plt
import argparse
from array import array

import sys
sys.path.append('..')
from SBCcode.DataHandling import ReadBinary as rb

from ozTools import *

import ROOT
from ROOT import TFile, TTree
from ROOT import TCanvas, gApplication, gStyle, TLegend, TPaveText
from ROOT import TF1, TH1F, TGraph, TGraphErrors, TFitResult
from ROOT import TSpectrum
from ROOT import TMath, TFormula
from ROOT import gROOT

print("Austin's makeshift static variable processor")

# index of recorded waveforms in the event
iSipm = 0 # sipm data is first
iDis = 1 # discriminator data second

###### parsing ######
pars = argparse.ArgumentParser()
pars.add_argument("--first", type = int, required = False, default = 0, help = "what event to start in")
pars.add_argument("--g", type = int, required = False, default = 0, help = "what pulse group to start on")
pars.add_argument("--N", type = int, required = False, default = 0, help = "N events to plot together, default is 0 for all events")
pars.add_argument("--run", type = str, required = False, default = "20210726_0")
pars.add_argument("--live", type = int, required = False, default = False, help = "do viewing with user input")
pars.add_argument("--view", type = int, required = False, default = False, help = "view final plots")
pars.add_argument("--outName", type = str, required = False)
args=pars.parse_args()
run = args.run
live = bool(args.live)
view = bool(args.view)

##### output and input file #####
savePath = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/"
raw_directory = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/rawdata"
outName = savePath + "outProcStatic1_" + str(run) + ".root"
if args.outName != None:
	outName = str(outPath + args.outName + str(run) + ".root")

hFile = gROOT.FindObject(outName)
if hFile:
	hFile.Close()
hFile = TFile(outName, 'RECREATE', 'bprocStaticVars1.py output file')
print("writing ROOT file to:", outName)

### defining TTree and all variables ###
tt = TTree("tt", "proc static variables, run:" + str(run))

height1 = np.empty((1), dtype = "float32")
tt.Branch("hbline_min", height1, "hbline_min/F")

blineVal = np.empty((1), dtype = "float32")
tt.Branch("blineV", blineVal, "blineV/F")

blineE = np.empty((1), dtype = "float32")
tt.Branch("blineE", blineE, "blineE/F")

qr1 = np.empty((1), dtype = "float32")
tt.Branch("qr300", qr1, "qr300/F")

qr2 = np.empty((1), dtype = "float32")
tt.Branch("qr500", qr2, "qr500/F")

qs1 = np.empty((1), dtype = "float32")
tt.Branch("qs50", qs1, "qs50/F")

f1q = np.empty((1), dtype = "float32")
tt.Branch("f1q", qF1, "f1q/F")

f1r = np.empty((1), dtype = "float32")
tt.Branch("f1r", f1r, "f1r/F")

f1f = np.empty((1), dtype = "float32")
tt.Branch("f1f", f1f, "f1f/F")

f1chi2 = np.empty((1), dtype = "float32")
tt.Branch("f1chi2", f1chi2, "f1chi2/F")

# a = np.empty(1, dtype = object)
# a[0] = TF1()
fpulse1 = TF1()
tt.Branch("fitPulse1", fpulse1)

for i in range(10):
	pulse1 = TF1("myPulse", "pol4")
	pulse1.SetParameter(0, i)
	tt.Fill()

# 
# tt.Fill()
hFile.Write()

hFile.Close()
sys.exit()



