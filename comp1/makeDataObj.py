
print("hi austin! Running makeDataObj")

## this is a shortened version of the analysis code written in other places.
# It loads in the data, does PE normalization, rate calculation and bkg subtraction

# it produces X output RC objects, whose histograms are useful for comparing to MC

# it is unclear how the mapping from rcObject to MC will be done
# user inputs must be automated
  # syke can handle all this in notebook

## ROOT stuff
import ROOT
from ROOT import TFile, TTree, TCanvas
from ROOT import TH1F, TH2F
import ctypes
import copy

ROOT.gStyle.SetCanvasDefW(600)
ROOT.gStyle.SetCanvasDefH(300)

# ROOT.gErrorIgnoreLevel = ROOT.kFatal
ROOT.gROOT.SetBatch(False)
ROOT.enableJSVis()
# ROOT.disableJSVis()

import os
import sys
import numpy as np

sys.path.append(os.environ.get("CF4HOME")+"/code1Analysis/")
sys.path.append(os.environ.get("CF4HOME")+"/code1Analysis/")

## Atown's stuff
from helpers.anaStructs import * # includes SORTING class/namespace
from helpers.getSPE import *
import helpers.graphBuilders as GB

# ------------- ------------- User Inputs  -------------------------- ##
## USER INPUTS (kinda bad ngl) ##
NCHAN = 3
VERBOSE = False

#QVAR_USE = "qRCFlat"
# QVAR_USE = "qRCFixed"
QVAR_USE = "qFlat"
# QVAR_USE = "qFixed"
PECUTVALUE = 0.5

## Paths to files, runLists
defPath = str(os.environ["CF4HOME"])+"/procData/" # where to get everything
bkgRunPath = "runLists/analysis1/bkgLists/"       # bkg run list files

## -------- physics and bkg runs -------- ##
runID_cs_bare_heights = [202203101341, 202203101306, 202203101150, 202203101534, 202203101506]
cs_bare_heights = [26.5, 28,	30.5,	32,	34.5]

runID_cs_lp_heights = [202204111600,	202204111454,	202204111423, 202204121321]
cs_lp_heights = [26, 28.5,	30.5, 32]

physicsList = runID_cs_bare_heights + runID_cs_lp_heights

# bkg runs
bkgList_0310 = np.loadtxt(defPath + bkgRunPath+ "date_0310-48V.txt", dtype=int).tolist()
bkgList_0411 = np.loadtxt(defPath + bkgRunPath+ "date_0411-48V.txt", dtype=int).tolist()

# final inputs into code
allBkgs = bkgList_0310 + bkgList_0411
allRuns = allBkgs + physicsList

## ------------- ------------- Getting RC objects ------------- ------------- ##
# Get all data RC objects, build needed SPE histograms ##
ROOT.gErrorIgnoreLevel = ROOT.kError
rcList = [] ## everything from here
fList = []

speList = []
chi2list = []

# loading data and SPE analysis
# everything is in rcList
for run in allRuns:
  f = ROOT.TFile(defPath + "analysis1/npv1_" + str(run) + ".root", "READ") #file opens
  fList.append(f) # NEED THIS to keep ttrees persistent (all reference in Python)

  tt = f.Get("tt")
  rcObj = RunClass(run, tt)
  rcObj.mFname = str(run)

  rcObj.assign_height(runID_cs_bare_heights, cs_bare_heights)
  rcObj.assign_height(runID_cs_lp_heights, cs_lp_heights)

  rcObj.computeSPE(width = 15, verbose=VERBOSE, varname=QVAR_USE, channels = [i for i in range(NCHAN)])
  rcList.append(rcObj)

  for chan in rcObj.mChannels:
    # print(chan.mSpe.chi2)
    speList.append(chan.mSpe.spe)
    chi2list.append(chan.mSpe.chi2)
ROOT.gErrorIgnoreLevel = ROOT.kUnset
print("Done Calculating SPE")

#### ------------- Sorting RC objects ------------- ####
rcl_Bkg_0310 = SORTING.pickListOfRC(bkgList_0310, rcList) # get all RC background objects
rcl_cs_bare_heights = SORTING.pickListOfRC(runID_cs_bare_heights, rcList)

rcl_Bkg_0411 = SORTING.pickListOfRC(bkgList_0411, rcList) # get all RC background objects
rcl_cs_lp_heights = SORTING.pickListOfRC(runID_cs_lp_heights, rcList)

######## ------------- Building Full Spectra -------------########
binLow = 0.5
binHigh = 100.5
Nbins = 100.
binWidth = (binHigh - binLow) / Nbins
hReference = TH1F("hReference", "hReference", int(Nbins), binLow, binHigh) # histogram that everything is built off
hReference.Sumw2(True)

# building all SPE spectra
for rc in rcList:
  rc.BuildSpectrum(hReference, QVAR_USE, PECUTVALUE, rate = True)

# building averaged bkg spectra
hNPE_bkg_0310 = hReference.Clone("hNPE_bkg_0310") # histogram
GB.build_bkgHistogramFast(rcl_Bkg_0310, QVAR_USE, hNPE_bkg_0310, PECUTVALUE) # this calls build_SPEFixedHistogram iteratively, does rate calculation...

hNPE_bkg_0411 = hReference.Clone("hNPE_bkg_0411") # histogram
GB.build_bkgHistogramFast(rcl_Bkg_0411, QVAR_USE, hNPE_bkg_0411, PECUTVALUE) # this calls build_SPEFixedHistogram iteratively, does rate calculation...

## ------------- Get and build the bkg corrected spectra ------------- ##
# print("\n Bare:")
# for rc in rcl_cs_bare_heights:
#   rc.FindClosestRC(rcl_Bkg_0310)
#   rc.BuildCorrSpectrum()
#   rc.setup_forComparison()
#   print(rc.mHeight)

# print("\n Filtered")
# for rc in rcl_cs_lp_heights:
#   rc.FindClosestRC(rcl_Bkg_0411)
#   rc.BuildCorrSpectrum()
#   rc.setup_forComparison()
#   print(rc.mHeight)