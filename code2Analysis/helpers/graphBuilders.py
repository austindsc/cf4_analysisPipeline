import numpy as np
import ROOT
from helpers.getSPE import *
# from helpers.anaStructs import RunClass

def myHistReset(hist):
  hist.Reset()
  hist.SetTitle(hist.GetName())
  return

# no type hinting for 'rc' beacuse it causes circular dependency
def build_SPEFixed_histogram(rc, varName: str, histRef: ROOT.TH1, peCut = 0.5):

  for ev in rc.tt:
    time = getattr(ev, "evTime")
    # print(type(time))
    qVar = getattr(ev, varName)
    chanPE = np.array([qVar[0]/rc.GetSPE(0), qVar[1]/rc.GetSPE(1), qVar[2]/rc.GetSPE(2)])
    evNPE = np.dot(chanPE, chanPE > peCut) #only select ones greater than peCur (second argument is array of 1s and 0s)

    if evNPE > peCut:
      histRef.Fill(evNPE)
  
  print("done filling histogram:", histRef.GetName())
  return

# call build_SPEFixed.., adds them all together.
def build_bkgHistogram(rcList:list, varName: str, histRef: ROOT.TH1, peCut = 0.5):
  
  Nentries = 0
  totalTime = 0.

  for rc in rcList:
    # rc.print()
    htemp = histRef.Clone(histRef.GetName()+str(rc.id))
    htemp.Reset()
    build_SPEFixed_histogram(rc = rc, varName=varName, histRef = htemp, peCut = peCut)
    # htemp.Scale(1/rc.GetRunTime())
    Nentries += htemp.GetEntries()
    totalTime += rc.GetRunTime()

    histRef.Add(histRef, htemp, 1., 1.)
    # build PE histogram from every backgrounds run. then use total ?

  ## do rate scaling at end
  histRef.Scale(1./totalTime)
  print("done building bkg histogram")
  return