
import os
import numpy as np
from array import array
import copy
import time

import ROOT
from ROOT import TFile, TTree, TCanvas
from ROOT import TH1F

# Atown's stuff
from helpers.anaStructs import *
from helpers.getSPE import *
from helpers.graphBuilders import *

## USER INPUTS
NCHAN = 3
VERBOSE = False
QVAR_USE = "qRCFlat"
PECUTVALUE = 0.5

defPath = str(os.environ["CF4HOME"])+"/procData/"
runPath = "runLists/analysis1/bkgLists/" + "date_0310-48V.txt"

## Actual function
if __name__ == '__main__':

  runList = np.loadtxt(defPath + runPath, dtype=str)

  rcList = []
  fList = []

  # loading data and SPE analysis
  for run in runList:
    f = ROOT.TFile(defPath + "analysis1/npv1_" + run + ".root", "READ") #file opens
    fList.append(f) # NEED THIS to keep ttrees persistent (all reference in Python)

    tt = f.Get("tt")
    rcObj = RunClass(run, tt)

    ## Calculating SPE values
    for ich in range(NCHAN):
      tempSpe = getSPEVval2(rcObj.tt, width = 15, verbose = VERBOSE, varName = QVAR_USE, channel = ich, runName = run)
      rcObj.getChan(ich).addSpe(tempSpe) 
    
    rcList.append(rcObj)


  # ## Got all files ##
  # pickRCevent = lambda runID: [rc for rc in rcList if rc.id == runID][0]
  # rc_cs_bkg_bare = pickRCevent(runID_cs_bkg_lp) #bad way to do search
  # rc_cs_bkg_bare2 = pickRCevent(runID_cs_bkg_bare2)
  # rc_cs_mid_bare = pickRCevent(runID_cs_mid_lp)
  
  rcList[0].tt.Show(0)

  # # bline = np.empty((3), dtype = "float32")
  # # rc_cs_bkg_bare.tt.SetBranchAddress("blineV", bline)
  # # for i in range(0, 10):
  # #   rc_cs_bkg_bare.tt.GetEntry(i)
  # #   print(bline)

  start = time.time()

  #### Building large histogram ####
  hNPE_bkg_0310 = TH1F("hNPE_bkg_0310", "hNPE_bkg_0310", 500, 0, 100)

  for rc in rcList:
    rc.print()
    htemp = hNPE_bkg_0310.Clone("htemp_NPE"+rc.id)
    htemp.Reset()
    build_SPEFixed_histogram(rc = rc, varName=QVAR_USE, histRef = htemp, peCut = PECUTVALUE)
    htemp.Scale(1/rc.GetRunTime())

    hNPE_bkg_0310.Add(hNPE_bkg_0310, htemp, 1., 1.)
    # build SPE histogram from every backgrounds run. then use total ?

    
  end = time.time()
  print ("Time elapsed:", end - start)
  
  # ########### Sanity Checks ###########
  # ## peak time 


  # ########### Building full PE histograms ###########

  # ########### Plotting ###########
  c1 = TCanvas()

  hNPE_bkg_0310.Draw("hist")

  c1.SetLogy()
  c1.Update()


  # # ROOT.gApplication.Run()
  # # pass