print("Running makeDataObj_co")

## this is a shortened version of the analysis code written in other places.
# It loads in the data, does PE normalization, rate calculation. leaves bkg subtraction for notebook (to select how its done)

# it produces X output RC objects, whose histograms are used for comparing to MC

# mapping from rc to mc is handled in notebook

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
import json5

sys.path.append(os.environ.get("CF4HOME")+"/code1Analysis/")
sys.path.append(os.environ.get("CF4HOME")+"/code1Analysis/")
THISPATH = os.environ.get("CF4HOME")+"/comp1/"

## Atown's stuff
from helpers.anaStructs import * # includes SORTING class/namespace
from helpers.getSPE import *
import helpers.graphBuilders as GB


# ------------- ------------- User Inputs  -------------------------- ##
with open(THISPATH+'config.json', 'r') as file:
  config = json5.load(file)

############## USER INPUTS ##############
NCHAN = config.get("data_analysis").get("nchan")
VERBOSE = config.get("data_analysis").get("verbose")

QVAR_USE = config.get("data_analysis").get("qvar")
PECUTVALUE = config.get("data_analysis").get("peCut")
SPEWIDTH = config.get("data_analysis").get("speWidth")

BINLOW = config.get("global_analysis").get("binLow")
BIGHIGH = config.get("global_analysis").get("binHigh")
NBINS = config.get("global_analysis").get("nBins")
binWidth = (BIGHIGH - BINLOW) / NBINS

## Paths to files, runLists
defPath = str(os.environ["CF4HOME"])+ config.get("data_analysis").get("defData_Path") # where to get everything
bkgRunPath = config.get("data_analysis").get("bkgLists_Path")       # bkg run list files

## -------- physics and bkg runs -------- ##
runID_co_bare_heights = config.get("cobaltBareData").get("runIDs")
co_bare_heights = config.get("cobaltBareData").get("heights")

runID_co_lp_heights = config.get("cobaltLPData").get("runIDs")
co_lp_heights = config.get("cobaltLPData").get("heights")

physicsList = runID_co_bare_heights + runID_co_lp_heights

# bkg runs
co_Bare_bkgFile = config.get("cobaltBareData").get("bkgList")
co_LP_bkgFile = config.get("cobaltLPData").get("bkgList")

# print(defPath, bkgRunPath, co_Bare_bkgFile)
bkgList_0303 = np.loadtxt(defPath + bkgRunPath + co_Bare_bkgFile, dtype=int).tolist()
bkgList_0413 = np.loadtxt(defPath + bkgRunPath + co_LP_bkgFile, dtype=int).tolist()

############## DONE USER INPUTS ##############

# final inputs into code
allBkgs = bkgList_0303 + bkgList_0413
allRuns = allBkgs + physicsList

## ------------- ------------- Getting RC objects ------------- ------------- ##
# Get all data RC objects, build needed SPE histograms ##
ROOT.gErrorIgnoreLevel = ROOT.kError
rcList = [] ## everything from here
fList = []

speListCo = []
chi2ListCo = []

# loading data and SPE analysis
# everything is in rcList
for run in allRuns:
  f = ROOT.TFile(defPath + "analysis1/npv1_" + str(run) + ".root", "READ") #file opens
  fList.append(f) # NEED THIS to keep ttrees persistent (all reference in Python)

  tt = f.Get("tt")
  rcObj = RunClass(run, tt)
  rcObj.mFname = str(run)

  rcObj.assign_height(runID_co_bare_heights, co_bare_heights)
  rcObj.assign_height(runID_co_lp_heights, co_lp_heights)

  rcObj.computeSPE(width = SPEWIDTH, verbose=VERBOSE, varname=QVAR_USE, channels = [i for i in range(NCHAN)])
  rcList.append(rcObj)

  for chan in rcObj.mChannels:
    # print(chan.mSpe.chi2)
    speListCo.append(chan.mSpe.spe)
    chi2ListCo.append(chan.mSpe.chi2)
ROOT.gErrorIgnoreLevel = ROOT.kUnset
print("Done Calculating SPE")

#### ------------- Sorting RC objects ------------- ####
rcl_Bkg_0303 = SORTING.pickListOfRC(bkgList_0303, rcList) # get all RC background objects
rcl_co_bare_heights = SORTING.pickListOfRC(runID_co_bare_heights, rcList)

rcl_Bkg_0413 = SORTING.pickListOfRC(bkgList_0413, rcList) # get all RC background objects
rcl_co_lp_heights = SORTING.pickListOfRC(runID_co_lp_heights, rcList)

######## ------------- Building Full Spectra -------------########
hReference = TH1F("hReference", "hReference", int(NBINS), BINLOW, BIGHIGH) # histogram that everything is built off
hReference.Sumw2(True)

# building all SPE spectra
for rc in rcList:
  rc.BuildSpectrum(hReference, QVAR_USE, PECUTVALUE, rate = True)

# building averaged bkg spectra
hNPE_bkg_0303 = hReference.Clone("hNPE_bkg_0303") # histogram
GB.build_bkgHistogramFast(rcl_Bkg_0303, QVAR_USE, hNPE_bkg_0303, PECUTVALUE) # this calls build_SPEFixedHistogram iteratively, does rate calculation...

hNPE_bkg_0413 = hReference.Clone("hNPE_bkg_0413") # histogram
GB.build_bkgHistogramFast(rcl_Bkg_0413, QVAR_USE, hNPE_bkg_0413, PECUTVALUE) # this calls build_SPEFixedHistogram iteratively, does rate calculation...


RAN_CO = True

jjjjjj = 69.420
## ------------- Get and build the bkg corrected spectra ------------- ##
# print("\n Bare:")
# for rc in rcl_co_bare_heights:
#   rc.FindClosestRC(rcl_Bkg_0303)
#   rc.BuildCorrSpectrum()
#   rc.setup_forComparison()
#   print(rc.mHeight)

# print("\n Filtered")
# for rc in rcl_co_lp_heights:
#   rc.FindClosestRC(rcl_Bkg_0413)
#   rc.BuildCorrSpectrum()
#   rc.setup_forComparison()
#   print(rc.mHeight)