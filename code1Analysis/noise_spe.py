## Python Stuff
import os
import sys
import numpy as np
from array import array
import copy
import time
import matplotlib.pyplot as plt

## ROOT stuff
import ROOT
from ROOT import TFile, TTree, TCanvas
from ROOT import TH1F, TH2F
import ctypes
import copy

ROOT.gStyle.SetCanvasDefW(600)
ROOT.gStyle.SetCanvasDefH(500)

# ROOT.gErrorIgnoreLevel = ROOT.kFatal
ROOT.gROOT.SetBatch(False)
# ROOT.enableJSVis()
# ROOT.disableJSVis()

sys.path.append(os.environ.get("CF4HOME")+"/code1Analysis/")
print(sys.path)

## Atown's stuff
from helpers.anaStructs import * # includes SORTING class/namespace
from helpers.getSPE import *
import helpers.graphBuilders as GB

## ------------- User Inputs ------------- ##
## USER INPUTS (kinda bad ngl) ##
NCHAN = 3
VERBOSE = True
ASK_MATPLOT = True
DO2D = False

# QVAR_USE = "qRCFlat"
# QVAR_USE = "qRCFixed"
QVAR_USE = "qFlat"
# QVAR_USE = "qFixed"
PECUTVALUE = 0.0

errorV = 0.05

## Paths to files, runLists
defPath = str(os.environ["CF4HOME"])+"/procData/" # where to get everything

## -------- physics and bkg runs -------- ##
bkgVoltagePath = "runLists/analysis1/"       # bkg run list files
bkgVoltageRuns = np.loadtxt(defPath + bkgVoltagePath + "bkg_voltage_runs.txt", dtype=int).tolist()
bkgVoltages = np.loadtxt(defPath + bkgVoltagePath + "bkg_voltage_list.txt", dtype=int).tolist()

# final inputs into code
allRuns = bkgVoltageRuns
allRuns

## Get all RC objects ##
rcList = [] ## everything from here
fList = []

speList = []
chi2list = []

# loading data and SPE analysis
# everything is in rcList
for run in allRuns[:]:
  f = ROOT.TFile(defPath + "analysis1/npv1_" + str(run) + ".root", "READ") #file opens
  fList.append(f) # NEED THIS to keep ttrees persistent (all reference in Python)

  tt = f.Get("tt")
  rcObj = RunClass(run, tt)
  rcObj.assign_voltage(allRuns, bkgVoltages)

  rcObj.computeSPE(width = 15, verbose=VERBOSE, varname=QVAR_USE, channels = [i for i in range(NCHAN)], nbins = 200, specialCase = rcObj.id)
  print(type(rcObj.id))
  rcObj.print()
  rcList.append(rcObj)

  # for chan in rcObj.mChannels:
  #   # print(chan.mSpe.chi2)
  #   speList.append(chan.mSpe.spe)
  #   chi2list.append(chan.mSpe.chi2)
if VERBOSE:
  input("say somthing to move past bkg SPE plots")


# check for outliers
ySPE = np.array([[rc.GetSPE(ich) for ich in range(NCHAN)] for rc in rcList], dtype=np.float64)
yeSPE = np.array([[rc.GetChan(ich).mSpe.espe for ich in range(NCHAN)] for rc in rcList], dtype=np.float64)
xVoltage = np.array([[rc.mVoltage for ich in range(NCHAN)] for rc in rcList], dtype=np.float64)
xeV_long = np.array([[errorV for ich in range(NCHAN)] for rc in rcList], dtype=np.float64)
print(xeV_long)

# line(xVoltage[:,0], slopeSPE, vbrList[0][0])
totalSPEGraph = ROOT.TGraphErrors(len(xeV_long.flatten()), xVoltage.flatten(), ySPE.flatten(), xeV_long.flatten(), yeSPE.flatten())
# fl_total = ROOT.TF1("fl_total", "[0]*(x - [1])")

# totalSPEGraph.Fit("fl_total", "", "", 30., 51.5)
# vGlobal = fl_total.GetParameters()[1]
# vEglobal = fl_total.GetParErrors()[1]

c1 = TCanvas("cwhat", "cwhat", 600, 400)
totalSPEGraph.Draw("AE*")

c1.Update()
c1.Draw()