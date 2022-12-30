import numpy as np
import ROOT

def build_SPEFixed_histogram(rc: RunClass, varName: str, histRef: ROOT.TH1):

  for ev in rc.tt:
    qVar = getattr(ev, varName)
    chanPE = np.array([qVar[0]/rc.GetSPE(0), qVar[1]/rc.GetSPE(1), qVar[2]/rc.GetSPE(2)])
    evNPE = np.dot(chanPE, chanPE > PECUTVALUE)

    if evNPE > PECUTVALUE:
      histRef.Fill(evNPE)
  
  print("done filling histogram:", histRef.GetName())
  return