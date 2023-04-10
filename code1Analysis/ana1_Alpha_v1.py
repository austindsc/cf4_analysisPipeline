## Python Stuff
import os
import numpy as np
from array import array
import copy
import time

import matplotlib.pyplot as plt

## ROOT stuff
import ROOT
from ROOT import TFile, TTree, TCanvas
from ROOT import TH1F, TH2F
# ROOT.gStyle.SetOptStat(0) # this breaks SPE verbosity

import ctypes
import copy

## Atown's stuff
from helpers.anaStructs import * # includes SORTING class/namespace
from helpers.getSPE import *
import helpers.graphBuilders as GB

## ------------- User Inputs ------------- ##
## USER INPUTS (kinda bad ngl) ##
NCHAN = 3
VERBOSE = True
ASK_MATPLOT = False
DO_2D = False
# QVAR_USE = "qRCFlat"
# QVAR_USE = "qRCFixed"
QVAR_USE = "qFlat"
# QVAR_USE = "qFixed"
PECUTVALUE = 0.5

## Paths to files, runLists
defPath = str(os.environ["CF4HOME"])+"/procData/" # where to get everything
rootPath = defPath+"/analysis1/"
# bkgRunPath = "runLists/analysis1/bkgLists/"       # bkg run list files

# For now just analyzing one run each. will look at temperature/pressure dependence later

# physics runs
# runID_LP = [202205201135, 202205201136, 202205201137]
# voltagesLP = [47., 48., 49]
# 47, 48 49 V
runID_lp = [202206231023]

# runID_bare = [202206211709, 202206220945]
# high P, low T
runID_bare = [202205201136]

physicsList = runID_bare + runID_lp

# # bkg runs
# bkgList_0310 = np.loadtxt(defPath + bkgRunPath+ "date_0310-48V.txt", dtype=int).tolist()
# bkgList_0411 = np.loadtxt(defPath + bkgRunPath+ "date_0411-48V.txt", dtype=int).tolist()

# final inputs into code
# allBkgs = bkgList_0310 + bkgList_0411
allRuns = physicsList

## ------------- Done User input ------------- ##

## Actual function
if __name__ == '__main__':

  rcList = [] ## everything from here
  fList = []

  # loading data and SPE analysis
  # everything is in rcList
  for run in allRuns:
    f = ROOT.TFile(rootPath + "npv1_" + str(run) + ".root", "READ") #file opens
    fList.append(f) # NEED THIS to keep ttrees persistent (all reference in Python)

    tt = f.Get("tt")
    rcObj = RunClass(run, tt)

    rcObj.computeSPE(width = 15, verbose=VERBOSE, varname=QVAR_USE, channels = [i for i in range(NCHAN)])
    rcList.append(rcObj)

  if VERBOSE:
    input("say somthing to move past bkg SPE plots")
  
  # printing information
  # rcList[0].tt.Show(0)
  for rc in rcList:
    rc.print()

  #### Sorting RC objects ####
  rcBare = SORTING.pickSingleRC(runID_bare[0], rcList)
  rcLp = SORTING.pickSingleRC(runID_lp[0], rcList)
  
  ########### Sanity Checks ###########
  ## peak time 
  # BORING

  ######## Building Full Spectra ########
  binHigh = 80.
  Nbins = 80.
  binWidth = (binHigh - 0.) / Nbins
  hReference = TH1F("hReference", "hReference", int(Nbins), 0., binHigh) # histogram that everything is built off
  hReference.Sumw2(True)

  # building all SPE spectra
  for rc in rcList:
    rc.BuildSpectrum(hReference, QVAR_USE, PECUTVALUE, rate = True)

  ## Time between events for one run
  hTimeBetweenEvents = TH1F("hTimeBetweenEvents", "hTimeBetweenEvents", 100, 0., 0.2)
  for rc in rcList:
    rc.BuildTTNE(hTimeBetweenEvents)

  ########### Plotting ###########
  # if VERBOSE == False:
  #   ROOT.gStyle.SetOptStat(0)

  ### 2d histograms ###
  
  if DO_2D:
    h2sharing = TH2F("h2sharing", "h2sharing", int(80), 0., binHigh, int(40), 0., 1.0)
    for rc in rcList:
      hthis = h2sharing.Clone("h2sharing_" + str(rc.id))
      hthis.SetTitle("h2sharing_" + str(rc.GetHeight()))
      GB.build_2d_foccupancyHistogram(rc, hthis, qvar = QVAR_USE)
      rc.mh2d = hthis

  cduh = []
  if DO_2D:
    for rc in rcList:
      c = TCanvas()
      cduh.append(c)
      c.cd()
      rc.mh2d.Draw("colz")
      c.SetLogz()
      c.Update()

    input("say somthing to continue past 2d plots")

  ## plotting all bkg spectra

  c0 = TCanvas()
  c0.Divide(1, 2, 0.01, 0.01, 0)

  
  hLp = rcLp.GetSpectrum()
  hLp.SetLineColor(2)
  hLp.Draw("hist E")

  rcBare.GetSpectrum().Draw("hist E same")

  c0.Update()
