## Python Stuff
import os
import numpy as np
from array import array
import copy
import time

## ROOT stuff
import ROOT
from ROOT import TFile, TTree, TCanvas
from ROOT import TH1F

import matplotlib.pyplot as plt

## Atown's stuff
from helpers.anaStructs import * # includes SORTING class/namespace
from helpers.getSPE import *
import helpers.graphBuilders as GB

## ------------- User Inputs ------------- ##
## USER INPUTS (kinda bad ngl) ##
NCHAN = 3
VERBOSE = False
QVAR_USE = "qRCFlat"
PECUTVALUE = 0.5

print("now here")

## Paths to files, runLists
defPath = str(os.environ["CF4HOME"])+"/procData/" # where to get everything
bkgRunPath = "runLists/analysis1/bkgLists/"       # bkg run list files

# physics runs
runID_cs_bare_heights = [202203101341, 202203101306, 202203101150, 202203101534, 202203101506]
cs_bare_heights = [26.5, 28,	30.5,	32,	34.5]

runID_cs_LP_heights = [202204111600,	202204111454,	202204111423]
cs_LP_heights = [26, 28.5,	30.5]

physicsList = runID_cs_bare_heights + runID_cs_LP_heights

# bkg runs
bkgList_0310 = np.loadtxt(defPath + bkgRunPath+ "date_0310-48V.txt", dtype=int).tolist()
bkgList_0411 = np.loadtxt(defPath + bkgRunPath+ "date_0411-48V.txt", dtype=int).tolist()

# final inputs into code
allBkgs = bkgList_0310 + bkgList_0411
allRuns = allBkgs + physicsList

## ------------- Done User input ------------- ##

## Actual function
if __name__ == '__main__':

  rcList = [] ## everything from here
  fList = []

  # loading data and SPE analysis
  # everything is in rcList
  for run in allRuns:
    f = ROOT.TFile(defPath + "analysis1/npv1_" + str(run) + ".root", "READ") #file opens
    fList.append(f) # NEED THIS to keep ttrees persistent (all reference in Python)

    tt = f.Get("tt")
    rcObj = RunClass(run, tt)

    rcObj.assign_height(runID_cs_bare_heights, cs_bare_heights)
    rcObj.assign_height(runID_cs_LP_heights, cs_LP_heights)

    rcObj.computeSPE(width = 15, verbose=VERBOSE, varname=QVAR_USE, channels = [i for i in range(NCHAN)])
    rcList.append(rcObj)

  # printing information
  rcList[0].tt.Show(0)
  for rc in rcList:
    rc.print()

  #### Sorting RC objects ####
  rcl_Bkg_0310 = SORTING.pickListOfRC(bkgList_0310, rcList) # get all RC background objects
  rcl_cs_bare_heights = SORTING.pickListOfRC(runID_cs_bare_heights, rcList)

  rcl_Bkg_0411 = SORTING.pickListOfRC(bkgList_0411, rcList) # get all RC background objects
  rcl_cs_lp_heights = SORTING.pickListOfRC(runID_cs_LP_heights, rcList)

  ########### Sanity Checks ###########
  ## peak time 
  # BORING

  ######## Building Full Spectra ########
  hReference = TH1F("hReference", "hReference", 500, 0, 50) # histogram

  # building all SPE spectra
  for rc in rcList:
    rc.BuildSpectrum(hReference, QVAR_USE, PECUTVALUE, rate = True)

  # building averaged bkg spectra
  hNPE_bkg_0310 = hReference.Clone("hNPE_bkg_0310") # histogram
  GB.build_bkgHistogramFast(rcl_Bkg_0310, QVAR_USE, hNPE_bkg_0310, PECUTVALUE) # this calls build_SPEFixedHistogram iteratively, does rate calculation...

  hNPE_bkg_0411 = hReference.Clone("hNPE_bkg_0411") # histogram
  GB.build_bkgHistogramFast(rcl_Bkg_0411, QVAR_USE, hNPE_bkg_0411, PECUTVALUE) # this calls build_SPEFixedHistogram iteratively, does rate calculation...
  
  # for rc in rcl_cs_bare_heights:
  #   rc.BuildCorrSpectrum(hNPE_bkg_0310)
  
  # for rc in rcl_cs_lp_heights:
  #   rc.BuildCorrSpectrum(hNPE_bkg_0411)
  
  ########### Plotting ###########
  
  ## Raw Spectra
  # c1 = TCanvas()
  
  # hNPE_bkg_0310.Draw("hist E")
  # # hNPE_bkg_0310.GetXaxis().SetRangeUser(0, 20)
  # rateBkg_0310 = hNPE_bkg_0310.Integral()

  # rcl_cs_bare_heights[1].getSpectrum().SetLineColor(2)
  # rcl_cs_bare_heights[1].getSpectrum().Draw("hist E same")
  
  # c1.SetLogy()
  # c1.Update()

  ## plotting all bkg spectra
  c0 = TCanvas()
  c0.Divide(1, 2)
  c0.cd(1).SetLogy()

  ROOT.gStyle.SetOptStat(0)

  hReference.GetXaxis().SetRangeUser(0, 20)
  hReference.Draw("hist")

  for rc in rcl_Bkg_0310:
    h = rc.GetSpectrum()
    h.SetLineColor(rcl_Bkg_0310.index(rc)+1)
    h.SetLineWidth(3)

    h.Draw("hist same E")


  c0.cd(2)
  c0.cd(2).SetLogy()
  hReference.Draw("hist")

  for rc in rcl_Bkg_0411:
    h = rc.GetSpectrum()
    h.SetLineColor(rcl_Bkg_0411.index(rc)+1)
    h.SetLineWidth(2)

    h.Draw("hist same E")

  c0.Update()

  #### corrected figs ####
  # c2 = TCanvas()

  # hReference.Draw("hist")
  # hReference.GetXaxis().SetRangeUser(0, 20.)
  # leg = ROOT.TLegend()
  # ROOT.gStyle.SetOptStat(0)

  # for rc in rcl_cs_bare_heights:
  #   hCorr = rc.getCorrSpectrum()

  #   hCorr.SetLineColor(rcl_cs_bare_heights.index(rc)+1)

  #   hCorr.Draw("hist same E")
  #   leg.AddEntry(hCorr, str(rc.mHeight), "L")

  # leg.Draw("same")
  # c2.SetLogy()
  # c2.Update()


  ## Matplotlib plotting ##
  # xSourceHeight = np.array([rc.GetHeight() for rc in rcl_cs_bare_heights])
  # xSource_EOR = np.array([rc.GetEORtime() for rc in rcl_cs_bare_heights])
  # yFixedRate = np.array([rc.GetRate() for rc in rcl_cs_bare_heights])
  # yFixedRateError = np.array([rc.GetRateError() for rc in rcl_cs_bare_heights])

  # # heightLabels = np.array([rc.mHeight() for rc in rcl_cs_bare_heights])

  # y0311_rate = np.array([rc.GetRate() for rc in rcl_Bkg_0310])
  # y0311_rate_error = np.array([rc.GetRateError() for rc in rcl_Bkg_0310])
  # x0311_fake = np.linspace(start = min(xSourceHeight), stop=max(xSourceHeight), num = len(rcl_Bkg_0310))
  # x0311_EOR = np.array([rc.GetEORtime() for rc in rcl_Bkg_0310])

  # fig, ax = plt.subplots(1, 1, figsize = (8, 4))
  # ax.errorbar(x = xSourceHeight, y = yFixedRate, yerr = yFixedRateError)
  # ax.errorbar(x0311_fake, y0311_rate, y0311_rate_error, marker = '.')
  # fig.show()

  # fig2, ax2 = plt.subplots(1, 1, figsize = (8, 4))
  # ax2.set_xlabel('EOR time')
  # ax2.errorbar(x = xSource_EOR, y = yFixedRate, yerr = yFixedRateError, linestyle = 'none')
  # ax2.errorbar(x0311_EOR, y0311_rate, y0311_rate_error, marker = '.', linestyle = 'none')

  # for i, height in enumerate(xSourceHeight):
  #   ax2.annotate(height, (xSource_EOR[i], yFixedRate[i]))
  # fig2.show()

  # ROOT.gApplication.Run()
  # pass