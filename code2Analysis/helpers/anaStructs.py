import copy
# add TF1 objects in here somewhere?

#add voltage to RunClass object

## ---------------------------------------------------------------- ##
class RunClass:

  def __init__(self, run, ttree, nChan = 3):
    self.id = str(run) #run ID (from integer)
    self.tt = ttree
    self.label=""
    self.mChannels = [Channel(i) for i in range(nChan)]
    return
  
  def getChan(self, index):
    return self.mChannels[index]

  def GetSPE(self, index):
    return self.mChannels[index].mSpe.spe

  def print(self):
    print("\nrunID:", self.id)
    for chan in self.mChannels:
      chan.print()
    return

## ---------------------------------------------------------------- ##
class Channel:
  def __init__(self, chID):
    self.mChid = chID
    self.mSpe = SpeStuff() # empty
    return
  
  def addSpe(self, spe):
    self.mSpe = spe
    return
  
  def print(self):
    print(f"chID:{self.mChid}:", end = " ")
    self.mSpe.print()
    return

## ---------------------------------------------------------------- ##
class SpeStuff:

  def __init__(self, spe = -1, speE = -1, chi2 = 0):
    self.spe = spe
    self.espe = speE
    self.chi2 = chi2 # maximum value
    return

  def print(self):
    print(f"spe = {self.spe:0.2f} +/- {self.espe:0.2f}     maxChi2 = {self.chi2:0.2f}")
    return