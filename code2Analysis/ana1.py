
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

## crappy run organization
runID_cs_bkg_bare = "202203101245"
runID_cs_mid_bare = "202203101306"

runID_cs_bkg_lp = "202204111434"
runID_cs_mid_lp = "202204111454"

## Actual function
if __name__ == '__main__':

  defPath = str(os.environ["CF4HOME"])+"/procData/"
  runList = np.loadtxt(defPath+ "runLists/analysis1/" + "csList.txt", dtype=str)

  rcList = []
  fList = []

  # loading data and SPE analysis
  for run in runList:
    f = ROOT.TFile(defPath + "analysis1/npv1_" + run + ".root", "READ") #file opens
    fList.append(f) # NEED THIS to keep ttrees persistent (all reference in Python)

    tt = f.Get("tt")
    rcObj = RunClass(run, tt)

    for ich in range(NCHAN):
      tempSpe = getSPEVval(rcObj.tt, width = 15, verbose = VERBOSE, varName = QVAR_USE, channel = ich, runName = run)
      rcObj.getChan(ich).addSpe(tempSpe) 
    
    rcList.append(rcObj)

  ## Got all files ##
  pickRCevent = lambda runID: [rc for rc in rcList if rc.id == runID][0]
  rc_cs_bkg_bare = pickRCevent(runID_cs_bkg_lp) #bad way to do search
  rc_cs_mid_bare = pickRCevent(runID_cs_mid_lp)
  
  rc_cs_bkg_bare.tt.Show(0)

  # bline = np.empty((3), dtype = "float32")
  # rc_cs_bkg_bare.tt.SetBranchAddress("blineV", bline)
  # for i in range(0, 10):
  #   rc_cs_bkg_bare.tt.GetEntry(i)
  #   print(bline)

  # rc_cs_bkg_bare.print()
  for rc in rcList:
    rc.print()

##xx finish print function. check all working. then check all SPE are good. then sanity checks.

  ########### Sanity Checks ###########
  ## peak time 

  ## 

  ########### Building full PE histograms ###########

  hNPE_csBkg_bare = TH1F("hNPE_bkg_bare", "hNPE_bkg_bare", 100, 0, 100)
  hNPE_csMid_bare = hNPE_csBkg_bare.Clone("hNPE_csMid_bare")
  hNPE_csMid_bare.Reset()

  start = time.time()

  build_SPEFixed_histogram(rc = rc_cs_bkg_bare, varName=QVAR_USE, histRef = hNPE_csBkg_bare)
  build_SPEFixed_histogram(rc = rc_cs_mid_bare, varName=QVAR_USE, histRef = hNPE_csMid_bare)
  
  end = time.time()
  print ("Time elapsed:", end - start)
  

  ########### Plotting ###########
  c1 = TCanvas()

  hNPE_csBkg_bare.Scale(1./hNPE_csBkg_bare.GetEntries())
  hNPE_csBkg_bare.Draw("hist")

  hNPE_csMid_bare.Scale(1./hNPE_csMid_bare.GetEntries())
  hNPE_csMid_bare.SetLineColor(2)
  hNPE_csMid_bare.Draw("hist same")

  c1.SetLogy()
  c1.Update()


  # ROOT.gApplication.Run()
  # pass