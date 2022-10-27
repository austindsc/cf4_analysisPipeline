import enum
import numpy as np

def doesThisCopy(a):
    a*=-1
    # a[0] = 69
    return

def howAboutThis(a):
    for i, val in enumerate(a):
        a[i] = -val
    return

def BaselineRestoration(wfm: np.array, tau: float, width: float = 16.0):
    
    vShift = np.zeros(len(wfm)) # initial condition
    corrWfm = np.zeros(len(wfm))
    
    factor = np.exp(-width/tau)
    print("factor = ", factor)

    # recursive shift
    for i, input in enumerate(wfm):
      #inpute = Vi
        correction = input - vShift[i]*factor #Vi'
        if(i+1 < len(wfm)):
            vShift[i+1] = vShift[i]*factor - correction*width/tau # Si+1
        
        corrWfm[i] = correction
    
    return corrWfm


def myRCrestore(wfm: np.array, tau: float, width: float = 16.0, bline: float = 0):

    corrWfm = wfm - bline
    vAdd = np.zeros(len(wfm))
    
    for i, V in enumerate(wfm):
      deltaV = wfm[i+1]-wfm[i] # value above baseline

      rcChangeToVoltage = -deltaV * width/tau

      corrWfm[i+1] += rcChangeToVoltage

    return


def undoRC(wfm, width, tau): 
  fixedRC = np.copy(wfm)

  for i, vp in enumerate(wfm[0:-1]):
    fixedRC[i+1] = (wfm[i+1] - wfm[i]) * tau/width + wfm[i]
  
  return fixedRC