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

## Atown's stuff
from helpers.anaStructs import *
from helpers.getSPE import *
from helpers.graphBuilders import *

## ------------- User Inputs ------------- ##
## USER INPUTS (kinda bad ngl) ##
NCHAN = 3
VERBOSE = False
QVAR_USE = "qRCFlat"
PECUTVALUE = 0.5

## Paths to files, runLists
defPath = str(os.environ["CF4HOME"])+"/procData/" # where to get everything
bkgRunPath = "runLists/analysis1/bkgLists/"       # bkg run list files

# physics runs
runID_cs_mid_bare = 202203101306
runID_cs_mid_lp = 202204111454
physicsList = [runID_cs_mid_bare, runID_cs_mid_lp]

# bkg runs
bkgList_0310 = np.loadtxt(defPath + bkgRunPath+ "date_0310-48V.txt", dtype=int).tolist()
bkgList_0411 = np.loadtxt(defPath + bkgRunPath+ "date_0411-48V.txt", dtype=int).tolist()
# runList = np.loadtxt(defPath+ "runLists/analysis1/" + "csList_noBkg.txt", dtype=str)

# final inputs into code
allBkgs = bkgList_0310 + bkgList_0411
allRuns = allBkgs + physicsList

## ------------- Done User input ------------- ##

## Actual function
if __name__ == '__main__':

  rcList = []
  fList = []

  # loading data and SPE analysis
  # everything is in rcList
  for run in allRuns:
    f = ROOT.TFile(defPath + "analysis1/npv1_" + str(run) + ".root", "READ") #file opens
    fList.append(f) # NEED THIS to keep ttrees persistent (all reference in Python)

    tt = f.Get("tt")
    rcObj = RunClass(run, tt)

    for ich in range(NCHAN):
      tempSpe = getSPEVval2(rcObj.tt, width = 15, verbose = VERBOSE, varName = QVAR_USE, channel = ich, runName = run)
      rcObj.getChan(ich).addSpe(tempSpe) 
    
    rcList.append(rcObj)

  # printing information
  rcList[0].tt.Show(0)
  for rc in rcList:
    rc.print()
  
  ########### Sanity Checks ###########
  ## peak time 
  # BORING


  ######## Building Full Spectra ########
  hReference = TH1F("hReference", "hReference", 500, 0, 100) # histogram

  ## 0310 data - Bare SiPM ##
  rcBkg_0310 = pickListOfRC(bkgList_0310, rcList) # get all RC background objects
  hNPE_bkg_0310 = TH1F("hNPE_bkg_0310", "hNPE_bkg_0310", 500, 0, 100) # histogram
  build_bkgHistogram(rcBkg_0310, QVAR_USE, hNPE_bkg_0310, PECUTVALUE) # this calls build_SPEFixedHistogram iteratively, does rate calculation...
  
  rc_cs_mid_bare = pickSingleRC(runID_cs_mid_bare, rcList) # physics run
  hNPE_cs_mid_bare = hNPE_bkg_0310.Clone("hNPE_cs_mid_bare")
  hNPE_cs_mid_bare.Reset()
  build_SPEFixed_histogram(rc_cs_mid_bare, QVAR_USE, hNPE_cs_mid_bare, PECUTVALUE)
  hNPE_cs_mid_bare.Scale(1./rc_cs_mid_bare.GetRunTime()) # rate scaling

  # bkg subtracted
  hNPE_cs_mid_bare_corrected = hNPE_cs_mid_bare.Clone("hNPE_cs_mid_bare_corrected")
  hNPE_cs_mid_bare_corrected.Add(hNPE_bkg_0310, -1.) # doing subtraction

  # test = PhysicsInstance()
  # test.mDate = "0310"
  # test.mBkgHist = hNPE_bkg_0310

  ## 0411 data - LP covered SiPM ##
  # repeat what is done above

  rcBkg_0411 = pickListOfRC(bkgList_0411, rcList)
  hNPE_bkg_0411 = hNPE_bkg_0310.Clone("hNPE_bkg_0411")
  myHistReset(hNPE_bkg_0411)
  build_bkgHistogram(rcBkg_0411, QVAR_USE, hNPE_bkg_0411, PECUTVALUE)

  rc_cs_mid_LP = pickSingleRC(runID_cs_mid_lp, rcList) # physics run
  hNPE_cs_mid_LP = hNPE_bkg_0310.Clone("hNPE_cs_mid_LP")
  myHistReset(hNPE_cs_mid_LP)
  build_SPEFixed_histogram(rc_cs_mid_LP, QVAR_USE, hNPE_cs_mid_LP, PECUTVALUE)
  hNPE_cs_mid_LP.Scale(1./rc_cs_mid_LP.GetRunTime()) 

  hNPE_cs_mid_LP_corrected = hNPE_cs_mid_LP.Clone("hNPE_cs_mid_LP_corrected")
  hNPE_cs_mid_LP_corrected.Add(hNPE_bkg_0411, -1.)

  ########### Plotting ###########
  
  ## Raw Spectra
  c1 = TCanvas("yo")

  c1.Divide(1, 2)

  c1.cd(1) # Bare data
  c1.cd(1).SetLogy()

  hNPE_bkg_0310.Draw("hist")

  hNPE_cs_mid_bare.SetLineColor(2)
  hNPE_cs_mid_bare.Draw("hist same")

  hNPE_cs_mid_bare_corrected.SetLineColor(3)
  hNPE_cs_mid_bare_corrected.SetLineWidth(4)
  hNPE_cs_mid_bare_corrected.Draw("hist same")

  c1.cd(2) # LP data
  c1.cd(2).SetLogy()

  hNPE_bkg_0411.Draw("hist")

  hNPE_cs_mid_LP.SetLineColor(2)
  hNPE_cs_mid_LP.Draw("hist same")

  hNPE_cs_mid_LP_corrected.SetLineColor(3)
  hNPE_cs_mid_LP_corrected.SetLineWidth(4)
  hNPE_cs_mid_LP_corrected.Draw("hist same")
  
  c1.Update()

  ## signal data

  c2 = TCanvas()

  # hNPE_cs_mid_bare_corrected.Rebin(2)
  # hNPE_cs_mid_LP_corrected.Rebin(2)


  hNPE_cs_mid_bare_corrected.SetLineColor(2)
  hNPE_cs_mid_bare_corrected.SetLineWidth(3)
  hNPE_cs_mid_bare_corrected.Draw("hist")

  hNPE_cs_mid_LP_corrected.SetLineColor(3)
  hNPE_cs_mid_LP_corrected.SetLineWidth(3)
  hNPE_cs_mid_LP_corrected.Draw("hist same")
  
  hNPE_cs_mid_bare_corrected.SetMaximum(1.1 * hNPE_cs_mid_LP_corrected.GetMaximum())
  hNPE_cs_mid_bare_corrected.SetMinimum(0.9 * 1e-4)
  hNPE_cs_mid_bare_corrected.GetXaxis().SetRangeUser(0, 20)

  hNPE_cs_mid_bare_corrected.SetTitle("Bkg subtracted - Bare vs LP - Liquid CF4 - Cs source")
  hNPE_cs_mid_bare_corrected.GetXaxis().SetTitle("PE")
  hNPE_cs_mid_bare_corrected.GetYaxis().SetTitle("Rate (Hz / some PE width)")

  c2.SetLogy()
  c2.Update()


  # ROOT.gApplication.Run()
  pass