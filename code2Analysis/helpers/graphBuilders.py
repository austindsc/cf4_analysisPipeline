import numpy as np
import ROOT
from helpers.anaStructs import *

def build_SPEFixed_histogram(rc: RunClass, varName: str, histRef: ROOT.TH1, peCut = 0.5):

  for ev in rc.tt:
    time = getattr(ev, "evTime")
    # print(type(time))
    qVar = getattr(ev, varName)
    chanPE = np.array([qVar[0]/rc.GetSPE(0), qVar[1]/rc.GetSPE(1), qVar[2]/rc.GetSPE(2)])
    evNPE = np.dot(chanPE, chanPE > peCut)

    if evNPE > peCut:
      histRef.Fill(evNPE)
  
  print("done filling histogram:", histRef.GetName())
  return