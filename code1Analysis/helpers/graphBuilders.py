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

# simpler, just adds up already calculated spectra of RC objects,
# scales by total number of entries. Averages rates (not adding them)
def build_bkgHistogramFast(rcList:list, varName: str, histRef: ROOT.TH1, peCut = 0.5):
  totalTime = 0.0
  for rc in rcList:
    if not hasattr(rc, 'mhSpectrum'):
      print("called build_bkgHistogramFast without rc object having mhSpectrum")
      exit(1)
    # weighted average for correct addition
    Noriginal = histRef.GetEntries()
    Nnew = rc.mhSpectrum.GetEntries()
    Ntot = Noriginal + Nnew
    histRef.Add(histRef, rc.mhSpectrum, Noriginal/Ntot, Nnew/Ntot)
    totalTime += rc.GetRunTime()

  ## do rate scaling at end
  print("done building fast bkg histogram for", histRef.GetName())
  print("total time:", totalTime, "\n")
  return

## Divide channel zero by 2
def build_2d_foccupancyHistogram(rc, th2, qvar = "qFixed", peCut = 0.5, norm = True):
  scaleCh0 = 2.
  if norm is False:
    scaleCh0 = 1.
  for ev in rc.tt:
    qVar = getattr(ev, qvar)
    chanPE = np.array([qVar[0]/rc.GetSPE(0), qVar[1]/rc.GetSPE(1), qVar[2]/rc.GetSPE(2)])
    evNPE = np.dot(chanPE, chanPE > 0.5) #only select ones greater than peCur (second argument is array of 1s and 0s)

    if evNPE > peCut:
      chanPEScaled = np.array([qVar[0]/scaleCh0/rc.GetSPE(0), qVar[1]/rc.GetSPE(1), qVar[2]/rc.GetSPE(2)])
      sharePE = max(chanPEScaled)/evNPE
      th2.Fill(evNPE, sharePE)
  
  print("built 2d sharing histogram for run", rc.id)
  return th2