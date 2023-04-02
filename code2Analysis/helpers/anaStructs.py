import copy
import numpy as np
import helpers.graphBuilders as GB
from helpers.getSPE import SpeStuff
import helpers.getSPE as GETSPE
# add TF1 objects in here somewhere?

## just some methods for sorting RC objects from global list
class SORTING:
  ## Got all files ##
  # pickSingleRC = lambda runID, rcList: [rc for rc in rcList if rc.id == runID][0]
  # get one rc event
  def pickSingleRC(runID, rcList):
    return [rc for rc in rcList if rc.id == runID][0]

  # pickListRC = lambda runID_list, rcList: [rc for rc in rcList if runID_list.count(rc.id) == 1]
  # select many rc events
  def pickListOfRC(runID_list, rcList):
    return  [rc for rc in rcList if runID_list.count(rc.id) == 1]


#add voltage to RunClass object

def getDate(runID):
	return str(runID)[4:4+4]

def getEORtime(runID):
  return str(runID)[-4:]

## ---------------------------------------------------------------- ##
class PhysicsInstance:
  def __init__(self):
    return

## ---------------------------------------------------------------- ##
class RunClass:

  def __init__(self, run, ttree, nChan = 3):
    self.id = run #run ID (from integer)
    self.tt = ttree
    self.label=""
    self.date = getDate(run)
    self.mChannels = [Channel(i) for i in range(nChan)]
    self.mHeight = 0 #0 means background
    return
  
  def GetChan(self, index):
    return self.mChannels[index]

  def GetSPE(self, index):
    return self.mChannels[index].mSpe.spe
  
  def GetRate(self):
    if hasattr(self, 'mRate'):
      return self.mRate
    else:
      self.mRate = self.tt.GetEntries()/self.GetRunTime()
      self.mRateError = np.sqrt(self.tt.GetEntries())/self.GetRunTime()
      return self.mRate
  
  def GetRateError(self):
    return self.mRateError
  
  def GetEORtime(self):
    # mins = getEORtime(self.runID)
    return int(getEORtime(self.id))

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
      
      # total time in 8ns bins = current time + number of rollovers*2^31
      self.mRunTime = 8*1e-9*(lastTime  + fullCount * rollOver) #units of seconds
    return self.mRunTime

  def print(self):
    print("\nrunID:", self.id)
    print("Nevents: %0.2f,   time: %0.2f min" % (self.tt.GetEntries(), self.GetRunTime() / 60))
    print("average Rate: %0.2f Hz" % self.GetRate())
    if self.mHeight != 0:
      print("height = %0.2f" % self.mHeight)

    for chan in self.mChannels:
      chan.print()
    return
  
  ## pick height from corresponding list of runs
  def assign_height(self, runIDList:list, runHeightList:list):
    if self.id in runIDList:
      self.mHeight = runHeightList[runIDList.index(self.id)]
    return
  
  def GetHeight(self):
    if hasattr(self, 'mHeight'):
      return self.mHeight
    print("no Height for rc", self.id)
    return None

  ## compute SPE for all channels
  def computeSPE(self, width = 15, verbose = True, varname = "qRCFlat", channels = [0, 1, 2]):
    for ich in channels:
      tempSPE = GETSPE.getSPEVval2(self.tt, width, verbose, varname, ich, str(self.id))
      self.GetChan(ich).addSpe(tempSPE)
    return

  ## Getting and building NPE spectra
  def getSpectrum(self):
    if hasattr(self, 'mhSpectrum'):
      return self.mhSpectrum
    else:
      print("rcObject" + self.id + " has no mhSpectrum")
      return None

  # if true, return PE vs rate spectrum (all channels summed)
  def BuildSpectrum(self, histClone, varname = "qRCFlat", peCut = 0.5, rate = True):
    self.mhSpectrum = histClone.Clone("hNPE_"+str(self.id))
    GB.build_SPEFixed_histogram(self, varname, self.mhSpectrum, peCut)
    if rate == True:
      self.mhSpectrum.Scale(1./self.GetRunTime())
    return self.mhSpectrum
  
  def GetSpectrum(self):
    return self.mhSpectrum
  
  def GetCorrSpectrum(self):
    if hasattr(self, 'mhCorrSpectrum'):
      return self.mhCorrSpectrum
    else:
      print("rcObject" + self.id + " has no mhCorrSpectrum")
      return None

  def BuildCorrSpectrum(self, histBkg):
    h = self.getSpectrum()
    hCorr = h.Clone(h.GetName()+"_corrected")
    hCorr.Add(histBkg, -1)
    self.mhCorrSpectrum = hCorr
    return hCorr
  
  ## time to next event histogram
  def BuildTTNE(self, histRef):
    evTime = np.array([0], dtype = "float32")
    self.tt.SetBranchAddress("evTime", evTime)
    self.tt.GetEntry(0)
    lastTime = evTime[0]
    for i in range(1, self.tt.GetEntries()):
      self.tt.GetEntry(i)
      newTime = evTime[0]
      deltaTime = (newTime-lastTime) * 8 * 1e-9
      if newTime < lastTime: # handle rollover
        deltaTime = 8 * 1e-9 * 2**31 - deltaTime
        lastTime = newTime
        # continue
      histRef.Fill(deltaTime)
      lastTime = newTime
    self.mhTTNE = histRef
    return self.mhTTNE
    


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

## SPE STUFF is in 'getSPE.py' ##