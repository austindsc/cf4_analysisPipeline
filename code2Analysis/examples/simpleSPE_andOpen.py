
import os
import numpy as np
from array import array
import ROOT
from ROOT import TFile, TTree, TCanvas
from ROOT import TH1F

import copy

from helpers.anaStructs import *
from helpers.getSPE import *

runID_cs_bkg_bare = "202203101245"
runID_cs_mid_bare = "202203101306"

runID_cs_bkg_lp = "202204111434"
runID_cs_mid_lp = "202204111454"


if __name__ == '__main__':

  defPath = str(os.environ["CF4HOME"])+"/procData/"
  runList = np.loadtxt(defPath+ "runLists/analysis1/" + "csList.txt", dtype=str)

  rcList = []
  fList = []
  for run in runList:
    f = ROOT.TFile(defPath + "analysis1/npv1_" + run + ".root", "READ") #file opens


    tt = f.Get("tt")
    rcObj = RunClass(run, tt)
    rcList.append(rcObj)
    fList.append(f) # NEED THIS
  
  print(rcList[2].tt.Show(0))

  for rc in rcList:
    rc.tt.Show(0)

  ## Got all files ##
  rc_cs_bkg_bare = [rc for rc in rcList if rc.id == runID_cs_bkg_bare][0]
  rc_cs_bkg_bare.tt.Show(0)

  # bline = np.empty((3), dtype = "float32")
  # rc_cs_bkg_bare.tt.SetBranchAddress("blineV", bline)
  # for i in range(0, 10):
  #   rc_cs_bkg_bare.tt.GetEntry(i)
  #   print(bline)

  # htemp = TH1F("htemp", "", 200, 1, 2000)
  # c1 = TCanvas()
  # rc_cs_bkg_bare.tt.Draw("qRCFlat>>htemp", "")

  # htemp.Draw("hist")
  # c1.Update()

  # c3 = TCanvas()
  # getSPEVval(rc_cs_bkg_bare.tt, width = 15, verbose = 1, varName = "qRCFlat", channel = 0)

  # # c4 = TCanvas()
  # getSPEVval(rc_cs_bkg_bare.tt, width = 15, verbose = 1, varName = "qRCFlat", channel = 1)

  for ch in range(3):
    getSPEVval(rc_cs_bkg_bare.tt, width = 15, verbose = 1, varName = "qRCFlat", channel = ch)

  # c1.Update()
  # c2.Update()
  # c1.SetLogy()
  # c1.Update()


  # ROOT.gApplication.Run()
  # pass