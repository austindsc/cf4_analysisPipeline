import numpy as np

def doesThisCopy(a):
    a*=-1
    # a[0] = 69
    return

def howAboutThis(a):
    for i, val in enumerate(a):
        a[i] = -val
    return

def BaselineRestoration(wfm: np.array, tau: float, t0: float = 0., width: float = 16.0):
    vShift = np.zeros(len(wfm)) # initial condition
    corrWfm = np.zeros(len(wfm))
    
    factor = np.exp(-width/tau)
    print("factor = ", factor)

    # recursive shift
    for i, input in enumerate(wfm):
        correction = input - vShift[i]*factor
        shift = vShift[i]*factor - correction*width/tau
        if(i+1 < len(wfm)):
            vShift[i+1] = shift
        
        corrWfm[i] = correction
    
    return corrWfm