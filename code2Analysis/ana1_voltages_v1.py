## Python Stuff
import os
import numpy as np
from array import array
import copy
import time


## ROOT stuff
import ROOT
from ROOT import TFile, TTree, TCanvas
from ROOT import TH1F, TH2F
# ROOT.gStyle.SetOptStat(0) # this breaks SPE verbosity

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

import ctypes
import copy

## Atown's stuff
from helpers.anaStructs import * # includes SORTING class/namespace
from helpers.getSPE import *
import helpers.graphBuilders as GB

## ------------- User Inputs ------------- ##
## USER INPUTS (kinda bad ngl) ##
NCHAN = 3
VERBOSE = False
# QVAR_USE = "qRCFlat"
# QVAR_USE = "qRCFixed"
QVAR_USE = "qFlat"
# QVAR_USE = "qFixed"
PECUTVALUE = 0.5

## Paths to files, runLists
defPath = str(os.environ["CF4HOME"])+"/procData/" # where to get everything
bkgRunPath = "runLists/analysis1/bkgLists/"       # bkg run list files

# bkg runs
runList_0414 = np.loadtxt(defPath + bkgRunPath+ "date_0414-manyV.txt", dtype=int).tolist()
voltageList = np.loadtxt(defPath + bkgRunPath+"0414_voltages.txt", dtype=int).tolist()

# removing last voltage (low stats)
runList_0414.pop()
voltageList.pop()

# final inputs into code
allRuns = runList_0414

## ------------- Done User input ------------- ##
  
def fvbr(x, p):
  return p[0]*(x[0] - p[1])

## Actual function
if __name__ == '__main__':

  rcList = [] ## everything from here
  fList = []

  # loading data and SPE analysis
  # everything is in rcList
  for run, vV in zip(allRuns, voltageList):
    f = ROOT.TFile(defPath + "analysis1/npv1_" + str(run) + ".root", "READ") #file opens
    fList.append(f) # NEED THIS to keep ttrees persistent (all reference in Python)

    tt = f.Get("tt")
    rcObj = RunClass(run, tt)

    rcObj.mVoltage = vV

    rcObj.computeSPE(width = 15, verbose=VERBOSE, varname=QVAR_USE, channels = [i for i in range(NCHAN)])
    rcList.append(rcObj)

  # printing information
  # rcList[0].tt.Show(0)
  for rc in rcList:
    rc.print()
  
  if VERBOSE:
    input("say somthing to move past bkg SPE/ROOT plots")
  else:
    [f.Close for f in fList]

  ######### Plotting #########

  ySPE = np.array([[rc.GetSPE(ich) for ich in range(NCHAN)] for rc in rcList])
  yeSPE = np.array([[rc.GetChan(ich).mSpe.espe for ich in range(NCHAN)] for rc in rcList])
  xV = np.array(voltageList, dtype=np.float64)
  errorV = 0.05
  xeV = np.repeat(errorV, len(xV))

  fitResults = []
  for ich in range(3):
    # ich = 0
    # Perform the linear fit
    ydata = np.array(ySPE[:,ich], dtype=np.float64)
    yerr = np.array(yeSPE[:,ich], dtype=np.float64)

    tg1 = ROOT.TGraphErrors(len(xV), xV, ydata, xeV, yerr)

    f1 = ROOT.TF1("f1_ch"+str(ich), fvbr, min(xV) - 5., max(xV) + 5., 2) # 2 parameters
    f1.SetParameters(0.1, 43)
    sFit = tg1.Fit(f1, "SR", "goff", 0., 50.)
    sFit.Print()
    fitResults.append(copy.deepcopy(sFit))

  plt.ion()
  fig, ax = plt.subplots(1, 1, figsize=(6, 4))

  
  # plt.show()