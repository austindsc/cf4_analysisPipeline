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
from ROOT import TGraphErrors, TLegend, RooFit, RooRealVar, RooDataSet, RooArgList, RooPlot, RooAddPdf

# ROOT.gStyle.SetOptStat(0) # this breaks SPE verbosity

# import matplotlib.pyplot as plt

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
  tgList = []
  for ich in range(3):
    # ich = 0
    # Perform the linear fit
    ydata = np.array(ySPE[:,ich], dtype=np.float64)
    yerr = np.array(yeSPE[:,ich], dtype=np.float64)

    tg1 = ROOT.TGraphErrors(len(xV), xV, ydata, xeV, yerr)
    tgList.append(tg1)

    f1 = ROOT.TF1("f1_ch"+str(ich), fvbr, min(xV) - 5., max(xV) + 5., 2) # 2 parameters
    f1.SetParameters(0.1, 43)
    sFit = tg1.Fit(f1, "SR", "goff", 0., 50.)
    sFit.Print()
    fitResults.append(copy.deepcopy(sFit))

  ## starting multifit
  y1 = np.array(ySPE[:,0], dtype=np.float64).tolist()
  ye1 = np.array(yeSPE[:,0], dtype=np.float64)
  tg1 = tgList[0]

  y2 = np.array(ySPE[:,1], dtype=np.float64).tolist()
  ye2 = np.array(yeSPE[:,1], dtype=np.float64)
  tg2 = tgList[1]

  y3 = np.array(ySPE[:,2], dtype=np.float64).tolist()
  ye3 = np.array(yeSPE[:,2], dtype=np.float64)
  tg3 = tgList[2]

  ### starting from 'combinedFit.py' example
  ipar1 = np.array([0, 1], dtype=np.int32) # shared slope and breakdown voltage 1
  ipar2 = np.array([0, 2], dtype=np.int32)
  ipar3 = np.array([0, 3], dtype=np.int32) # shared slope and breakdown voltage 2

  class GlobalChi2(object):
      def __init__(self, f1, f2, f3):
          self._f1 = f1
          self._f2 = f2
          self._f3 = f3

      def __call__(self, par):
          # parameter vector is first background (in common 1 and 2) and then is
          # signal (only in 2)

          # the zero-copy way to get a numpy array from a double *
          par_arr = np.frombuffer(par, dtype=np.float64, count=-1)

          p1 = par_arr[ipar1]
          p2 = par_arr[ipar2]
          p3 = par_arr[ipar3]

          return self._f1(p1) + self._f2(p2) + self._f3(p3)

  ## wtf ##
  fl1 = ROOT.TF1("fl1", "[0]*(x - [1])", 40, 60)
  fl2 = ROOT.TF1("fl2", "[0]*(x - [1])", 40, 60) # last is dummy
  fl3 = ROOT.TF1("fl3", "[0]*(x - [1])", 40, 60)

  wf1 = ROOT.Math.WrappedMultiTF1(fl1, 1)
  wf2 = ROOT.Math.WrappedMultiTF1(fl2, 1)
  wf3 = ROOT.Math.WrappedMultiTF1(fl3, 1)
  opt = ROOT.Fit.DataOptions()

  range1 = ROOT.Fit.DataRange()
  range1.SetRange(0, 100)
  data1 = ROOT.Fit.BinData(opt, range1)
  ROOT.Fit.FillData(data1, tg1)

  range2 = ROOT.Fit.DataRange()
  range2.SetRange(0, 100)
  data2 = ROOT.Fit.BinData(opt, range2)
  ROOT.Fit.FillData(data2, tg2)

  range3 = ROOT.Fit.DataRange()
  range3.SetRange(0, 100)
  data3 = ROOT.Fit.BinData(opt, range3)
  ROOT.Fit.FillData(data3, tg3)

  chi1 = ROOT.Fit.Chi2Function(data1, wf1)
  chi2 = ROOT.Fit.Chi2Function(data2, wf2)
  chi3 = ROOT.Fit.Chi2Function(data3, wf3)
  globalChi2 = GlobalChi2(chi1, chi2, chi3)

  fitter = ROOT.Fit.Fitter()

  # slope slope offset blah
  par0 = np.array([35, 43, 43, 43], dtype=np.float64)
  buffpar0 = np.frombuffer(par0, dtype=np.float64, count=-1)

  fitter.Config().SetParamsSettings(4, buffpar0)

  # fitter.Config().ParSettings(3).Fix()
  # fitter.Config().ParSettings(1).SetLimits(0, 50)
  # fitter.Config().ParSettings(2).SetLimits(0, 30)
  # fitter.Config().ParSettings(4).SetLimits(50, 100)
  # fitter.Config().ParSettings(5).SetLimits(0, 30)

  fitter.Config().MinimizerOptions().SetPrintLevel(0)
  fitter.Config().SetMinimizer("Minuit2", "Migrad")

  globalChi2Functor = ROOT.Math.Functor(globalChi2, 4)

  # # fit FCN function
  # # (specify optionally data size and flag to indicate that is a chi2 fit)
  fitter.FitFCN(globalChi2Functor, 0, data1.Size() + data2.Size() + data3.Size(), True)
  result = fitter.Result()
  result.Print(ROOT.std.cout)

  # c1 = ROOT.TCanvas("Simfit", "Simultaneous fit of two histograms", 10, 10, 700, 700)
  # c1.Divide(1, 2)
  # c1.cd(1)
  # ROOT.gStyle.SetOptFit(1111)

  # ftest.SetFitResult(result, ipar1)
  # ftest.SetRange(range1().first, range1().second)
  # ftest.SetLineColor(ROOT.kBlue)
  # ht.GetListOfFunctions().Add(ftest)
  # ht.Draw()

  # c1.cd(2)
  # ftest2.SetFitResult(result, ipar2)
  # ftest2.SetRange(range2().first, range2().second)
  # ftest2.SetLineColor(ROOT.kRed)
  # ht2.GetListOfFunctions().Add(ftest2)
  # ht2.Draw()

  # ht2.Draw()
  # c1.Update()