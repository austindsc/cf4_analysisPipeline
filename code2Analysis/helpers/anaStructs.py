import copy
import numpy as np
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

  def GetRunTime(self): # returns in seconds
    if not hasattr(self, 'mRunTime'):
      evTime = np.array([0], dtype = "float32")
      self.tt.SetBranchAddress("evTime", evTime)
      self.tt.GetEntry(0)
      totalTime = 0
      rollOver = 0

      fullCount = 2**31 - 1
      lastTime = int(evTime[0]) # got first time

      for i in range(1, self.tt.GetEntries()):
        self.tt.GetEntry(i)
        newTime = int(evTime[0])

        # check for rollover
        if newTime < lastTime:
          rollOver += 1
        lastTime = int(evTime[0])
      
      # total time in 16ns bins = current time + number of rollovers*2^31
      self.mRunTime = 16*1e-9*(lastTime  + fullCount * rollOver) #units of seconds
    return self.mRunTime

  def print(self):
    print("\nrunID:", self.id)
    print("Nevents: %0.2f,   time: %0.2f min" % (self.tt.GetEntries(), self.GetRunTime() / 60))
    print("average Rate: %0.2f Hz" % (self.tt.GetEntries()/self.GetRunTime()))
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