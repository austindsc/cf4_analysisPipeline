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


## USER INPUTS
NCHAN = 3
VERBOSE = False
QVAR_USE = "qRCFlat"
PECUTVALUE = 0.5

## crappy run organization
runID_cs_bkg_bare = "202203101333"
runID_cs_mid_bare = "202203101306"

runID_cs_bkg_lp = "202204111434"
runID_cs_mid_lp = "202204111520"


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

    # for ich in range(NCHAN):
    #   tempSpe = getSPEVval(rcObj.tt, width = 15, verbose = VERBOSE, varName = QVAR_USE, channel = ich, runName = run)
    #   rcObj.getChan(ich).addSpe(tempSpe) 
    
    rcList.append(rcObj)

  ## Got all files ##
  pickRCevent = lambda runID: [rc for rc in rcList if rc.id == runID][0]
  rc_cs_bkg_bare = pickRCevent(runID_cs_bkg_lp) #bad way to do search
  rc_cs_mid_bare = pickRCevent(runID_cs_mid_lp)
  
  rc_cs_mid_bare.tt.Show(0)

  counter = 0
  for ev in rc_cs_mid_bare.tt:
    # print(ev.evTime)
    if ev.evTime % 32 != 0:
      print(ev.evTime)
    # print(ev.evTime % 32, end = '\n')

    counter += 1
    if counter > 1000:
      break

  bevTime = np.array(1*[0], dtype = "float32")
  rc_cs_mid_bare.tt.SetBranchAddress("evTime", bevTime)
  rc_cs_mid_bare.tt.GetEntry(0)
  # rc_cs_mid_bare.tt.Print()
  # print(f"{bevTime[0]:0.10f}")

  rc_cs_bkg_bare.GetRunTime()
  rc_cs_bkg_bare.print()