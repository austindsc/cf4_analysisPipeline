"""
anaStructs.py - Data Analysis Structures for CF4 SBC Analysis

This module provides the core data structures for analyzing CF4 Super-heated Bubble Chamber data.
It handles ROOT files, event processing, SPE (Single Photo-Electron) analysis, and spectral analysis.

=== CLASSES OVERVIEW ===

1. SORTING (Static Class)
   - Utility methods for finding and filtering RunClass objects from lists
   - pickSingleRC(runID, rcList): Get single run by ID
   - pickListOfRC(runID_list, rcList): Get multiple runs by ID list  
   - findRCbyHeight(h, rcList): Find run by height parameter

2. RunClass (Main Analysis Class)
   - Represents a single data run with associated TTree and analysis results
   - Handles event timing, rate calculations, SPE analysis, and spectrum building
   
   Key Methods:
   - __init__(run, ttree, nChan=3, fname=""): Initialize with run ID and ROOT TTree
   - computeSPE(): Calculate Single Photo-Electron values for all channels
   - BuildSpectrum(): Create PE (Photo-Electron) vs rate spectrum
   - BuildCorrSpectrum(): Create background-subtracted spectrum
   - GetRate(), GetRunTime(): Calculate event rates and run duration
   - UprootConstructor(): Set up uproot interface for numpy-friendly analysis
   - GetSPEValues(), GetSPEArray(): Convert charge to SPE-reduced values
   - BuildCorrectEventTimes(): Handle ADC time rollover corrections
   - setup_forComparison(): Prepare data for matplotlib plotting

3. Channel (Per-Channel Data)
   - Represents individual PMT channel data
   - Stores SPE calibration information
   - addSpe(): Add SPE calibration data
   - print(): Display channel information

=== TYPICAL WORKFLOW ===
1. Create RunClass objects from ROOT files
2. Compute SPE calibrations with computeSPE()
3. Build event spectra with BuildSpectrum()
4. Perform background subtraction with BuildCorrSpectrum()
5. Use SORTING methods to organize and filter runs
6. Export to numpy arrays for further analysis

=== DEPENDENCIES ===
- ROOT/PyROOT for histogram and tree handling
- uproot for numpy-friendly ROOT file access
- helpers.graphBuilders for spectrum construction
- helpers.getSPE for SPE calibration algorithms
"""

import copy
import numpy as np
import sys
import helpers.graphBuilders as GB
from helpers.getSPE import SpeStuff
import helpers.getSPE as GETSPE
# add TF1 objects in here somewhere?

import uproot

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
  
  def findRCbyHeight(h, rcList):
    # FindRCbyHeight = lambda h, rclist: next((rc for rc in rclist if rc.GetHeight() == h), None)
    found = next((rc for rc in rcList if rc.GetHeight() == h), None)
    if found == None:
      print(f"ERROR, no rc object found for height {h}")
    return found


def getDate(runID):
	return str(runID)[4:4+4]

## ---------------------------------------------------------------- ##
class RunClass:

  def __init__(self, run, ttree, nChan = 3, fname = ""):
    self.id = run #run ID (from integer)
    self.tt = ttree
    self.label=""
    self.mDate = getDate(run)
    self.mChannels = [Channel(i) for i in range(nChan)]
    self.mHeight = 0 #0 means background
    self.mFname = fname
    self.mNEvents = int(ttree.GetEntries()) if ttree is not None else -1
    # self.mNEvents = int(69)
    return
  
  ## for averaging
  @classmethod
  def copy_rc_empty(cls, rcIn, name=""):
    newRc = cls(run = -1, ttree = None, nChan = 3, fname = name)
    newRc.mBin_edges = copy.deepcopy(rcIn.mBin_edges)
    newRc.mBinX = copy.deepcopy(rcIn.mBinX)
    newRc.mBin_rawCounts = np.zeros_like(rcIn.mBin_rawCounts)
    return newRc
  
  ## more averaging
  def AverageRCs(self, rcList):
    totalRunTime = 0
    for rc in rcList:
      self.mBin_rawCounts += rc.mBin_rawCounts
      totalRunTime += rc.GetRunTime()
    
    self.mBin_rate = self.mBin_rawCounts / totalRunTime
    self.mBin_rateError = np.sqrt(self.mBin_rawCounts) / totalRunTime
    return
  
  @classmethod
  def CombineRCs(cls, rcList, name="combined"):
    ## check that all have same binning
    for rc in rcList:
      if not np.array_equal(rc.mBin_edges, rcList[0].mBin_edges):
        raise ValueError("All RunClass objects must have the same binning to combine.")
    newRc = cls.copy_rc_empty(rcList[0], name=name)
    newRc.AverageRCs(rcList)
    return newRc
  
  def BackgroundSubtraction(self, rcbkg):
    ## existing members are: mBin_edges, mBinX, mBin_rawCounts, mBin_rate, mBin_rateError
    
    Nsignal = self.mBin_rawCounts - rcbkg.mBin_rate * self.GetRunTime()
    sigmaBkg = np.sqrt( (self.GetRunTime() * rcbkg.mBin_rateError)**2 + rcbkg.mBin_rate * self.GetRunTime() )
    error_signal = np.sqrt( Nsignal + sigmaBkg**2 )
  
    self.BkgCorrected_rawCounts = Nsignal
    self.BkgCorrected_errors = error_signal
    self.BkgCorrected_rate = Nsignal / self.GetRunTime()
    self.BkgCorrected_rateError = error_signal / self.GetRunTime()
    return
  
  ########## root stuff below (until np mentioned) ##########
  def setup_forComparison(self):
    # do a bunch of stuff to get some bins. moving to ROOT agnostic approach (but requires root spectrum)
    # for comparison to data
    if hasattr(self, "mhCorrSpectrum"):
      self.mBin_rate = np.array([self.mhCorrSpectrum.GetBinContent(i) for i in range(1, self.mhCorrSpectrum.GetNbinsX() + 1)])
      self.mBin_rawCounts = self.mBin_rate * self.GetRunTime()
      self.mBin_edges = np.array([self.mhCorrSpectrum.GetBinLowEdge(i) for i in range(1, self.mhCorrSpectrum.GetNbinsX() + 2)])
      self.mBinX = self.mBin_edges[:-1]
    elif hasattr(self, "mhSpectrum"):
      self.mBin_rate = np.array([self.mhSpectrum.GetBinContent(i) for i in range(1, self.mhSpectrum.GetNbinsX() + 1)])
      self.mBin_rawCounts = self.mBin_rate * self.GetRunTime()
      self.mBin_edges = np.array([self.mhSpectrum.GetBinLowEdge(i) for i in range(1, self.mhSpectrum.GetNbinsX() + 2)])
      self.mBinX = self.mBin_edges[:-1]
    return True
  
  def GetChan(self, index):
    return self.mChannels[index]

  def GetSPE(self, index):
    return self.mChannels[index].mSpe.spe
  
  def GetRate(self): # in Hz
    if hasattr(self, 'mRate'):
      return self.mRate
    else:
      self.mRate = self.mNEvents/self.GetRunTime()
      self.mRateError = np.sqrt(self.mNEvents)/self.GetRunTime()
      return self.mRate
  
  def GetRateError(self):
    return self.mRateError
  
  def GetEORtime(self):
    return str(self.id)[-4:]
  
  def GetEORmin(self):
    time = self.GetEORtime()
    return float(time[0:2])*60 + float(time[2:4])

  def GetRunTime(self): # returns in seconds
    if not hasattr(self, 'mRunTime'):
      evTime = np.array([0], dtype = "float32")
      self.tt.SetBranchAddress("evTime", evTime)
      self.tt.GetEntry(0)
      totalTime = 0
      rollOver = 0

      fullCount = 2**31 - 1
      lastTime = int(evTime[0]) # got first time

      for i in range(1, self.mNEvents):
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
    print("Nevents: %0.2f,   time: %0.2f min" % (self.mNEvents, self.GetRunTime() / 60))
    print("average Rate: %0.2f Hz" % self.GetRate())
    if self.mHeight != 0:
      print("height = %0.2f" % self.mHeight)

    for chan in self.mChannels:
      chan.print()
    return
  
  def FindClosestRC(self, rcl_bkg):
    t0 = self.GetEORmin()
    timeDifs = [abs(t0 - rc.GetEORmin()) for rc in rcl_bkg]
    indexClosestRC = timeDifs.index(min(timeDifs))
    self.mBkgRC = rcl_bkg[indexClosestRC]
    return self.mBkgRC
  
  ## pick height from corresponding list of runs
  def assign_height(self, runIDList:list, runHeightList:list):
    if self.id in runIDList:
      self.mHeight = runHeightList[runIDList.index(self.id)]
    return
  
  def assign_voltage(self, runIDList:list, voltageList:list):
    if self.id in runIDList:
      self.mVoltage = voltageList[runIDList.index(self.id)]
    return
  
  def GetHeight(self):
    if hasattr(self, 'mHeight'):
      return self.mHeight
    print("no Height for rc", self.id)
    return None

  ## compute SPE for all channels
  def computeSPE(self, width = 15, verbose = True, varname = "qFlat", channels = [0, 1, 2], nBins = 100, specialCase = 0):
    for ich in channels:
      tempSPE = GETSPE.getSPEVval2(self.tt, width, verbose, varname, ich, str(self.id), nBins = nBins, specialCase = specialCase)
      # tempSPE = GETSPE.getSPEVval(self.tt, width, verbose, varname, ich, str(self.id))
      self.GetChan(ich).addSpe(tempSPE)
    return

  ## Getting and building NPE spectra
  def GetSpectrum(self):
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
  
  def GetCorrSpectrum(self):
    if hasattr(self, 'mhCorrSpectrum'):
      return self.mhCorrSpectrum
    else:
      print("rcObject" + str(self.id) + " has no mhCorrSpectrum")
      return None

  ## pass it a bkg spectrum to subtract
  # call with nothing, it uses the pre-existing assigned spectrum from
  # the closest bkg RC object
  def BuildCorrSpectrum(self, histBkg = None):
    if histBkg is None:
      histBkg = self.mBkgRC.GetSpectrum()
    h = self.GetSpectrum()
    hCorr = h.Clone(h.GetName()+"_corrected")
    hCorr.Add(histBkg, -1)
    self.mhCorrSpectrum = hCorr

    self.mCorrBkgRate = histBkg.Integral() # store the rate bkg data
    self.Ncorr = hCorr.GetEntries()
    return hCorr
  
  ## time to next event histogram
  def BuildTTNE(self, histRef):
    self.mhTTNE = histRef.Clone(histRef.GetName() + "_" + str(self.id))
    evTime = np.array([0], dtype = "float32")
    self.tt.SetBranchAddress("evTime", evTime)
    self.tt.GetEntry(0)
    lastTime = evTime[0]
    for i in range(1, self.mNEvents):
      self.tt.GetEntry(i)
      newTime = evTime[0]
      deltaTime = (newTime-lastTime) * 8 * 1e-9
      if newTime < lastTime: # handle rollover
        deltaTime = 8 * 1e-9 * 2**31 - deltaTime
        lastTime = newTime
        # continue
      self.mhTTNE.Fill(deltaTime)
      lastTime = newTime
    return self.mhTTNE
  
  # def BuildTTNE_spe(self, histRef):
  #   self.m2dTTNE = histRef.Clone(histRef.GetName() + "_" + str(self.id))
  #   uprootTree = uproot.from_TTree(self.tt)
  #   pulseTimes = uprootTree["timeMin"]

###### methods below here are root agnostic ######
  def UprootConstructor(self, treeName = "tt"):
    fup = uproot.open(self.mFname)
    self.mttUp = fup["tt"]
    return

  ## convert tree charge values to SPE reduced values
  def GetSPEValues(self, qvar):
      self.npQPE_Arr = self.mttUp[qvar].array(library="np")/np.array([self.GetSPE(0), self.GetSPE(1), self.GetSPE(2)])
      return self.npQPE_Arr

  # this is the sum of SPE values for all channels
  def GetSPEArray(self, qvar="qRCFlat"):
      if not hasattr(self, "npQPE_Arr"):
          self.GetSPEValues(qvar)
      return np.sum(self.npQPE_Arr, axis=1)

  def BuildCorrectEventTimes(self): # using uproot, returns in ns
    if not hasattr(self, "mttUp"):
      print("ERROR: doesn't have uproot TTree")
      sys.exit(1)
    timeConvert = 8 # one clock's ADC value is 8 ns for event
    rollTime =  (2**31)* timeConvert # after hitting 2^32-1, it goes back to time + 2^31 (back to zero, except the last bit is flipped).
    
    evTime = np.array(self.mttUp["evTime"].array()) * timeConvert
    overIndexes = np.where(evTime[1:] < evTime[:-1])[0] + 1 # where is it larger than the previous entry
    mask = np.zeros_like(evTime)
    mask[overIndexes] = 1
    mask = np.cumsum(mask)
    self.mTrueEvTime = evTime + rollTime * mask
    return self.mTrueEvTime

  # need to call 'BuildCorrectEventTimes() first
  # this builds the pulse time within event and adds it to the true event time
  # returns timeMin with PE cut applied (in ns) and NPElist (for all events)
  def BuildCorrectEventPulseTimes(self, qvar, peCut = 0.5): # in ns
    # data is from my analysis code so times already in ns
    if not hasattr(self, "mttUp"):
      print("ERROR: doesn't have uproot TTree")
      sys.exit(1)
    pulseTimes = np.array(self.mttUp["timeMin"].array()) # time of waveform minimum (pulse peak)
    qList = np.array(self.mttUp[qvar].array())
    chanPEList = qList/np.array([self.GetSPE(0), self.GetSPE(1), self.GetSPE(2)])
    # chanPEList = np.sum(chanPEList, axis=1)
  
    # filter for PE > 0.5
    # peTimeMask = chanPEList > peCut
    # # print(peTimeMask[0:10])
    # timeMinmasked = time * peTimeMask
    # timeMinmasked[timeMinmasked == 0] = np.inf
    # timeMin_halfPE = np.min(timeMinmasked, axis=1)
    # timeMin_halfPE[timeMin_halfPE == np.inf] = 0 # fix previous step
    # self.mTimeMin_peCut = timeMin_halfPE
    
    return self.mTrueEvTime, pulseTimes, chanPEList
  
## methods adding since starting alpha analysis ##
  def setup_forComparison_np(self, Nbins = 100, lowRange=0.5, maxRange=100.5, qvar="qFlat", peCut=None):
    """
    ROOT-agnostic setup for comparison using numpy arrays directly.
    Creates same members as setup_forComparison() but bypasses mhSpectrum.
    
    Args:
        binArgs: Dictionary with binning arguments
                 {'nbins': int, 'low': float, 'high': float} or
                 {'edges': np.array} for custom bin edges
        qvar: Charge variable name (default "qFlat")
        peCut: PE threshold cut (default 0.5)
    """
    # Get PE array using existing method
    pe_array = self.GetSPEArray(qvar)
    
    # Apply PE cut
    pe_filtered = pe_array[pe_array >= peCut] if peCut is not None else pe_array

    # Create bin edges
    self.mBin_edges = np.linspace(lowRange, maxRange, Nbins + 1)

    # Create histogram - FIX: Use explicit range to exclude underflow/overflow
    self.mBin_rawCounts, _ = np.histogram(pe_filtered, bins=self.mBin_edges, 
                                         range=(self.mBin_edges[0], self.mBin_edges[-1]))
    
    # Calculate rates and bin centers
    self.mBin_rate = self.mBin_rawCounts / self.GetRunTime()
    self.mBin_rateError = np.sqrt(self.mBin_rawCounts) / self.GetRunTime()
    self.mBinX = self.mBin_edges[:-1] + np.diff(self.mBin_edges) / 2
    
    return True
  
  ## plotting helpers
  def mpl_hist(self, ax, rate = False, bkgCorrected = False, **kwargs):
    weights = self.mBin_rawCounts if rate == False else self.mBin_rate
    
    if bkgCorrected == True and hasattr(self, "BkgCorrected_rawCounts"):
      weights = self.BkgCorrected_rate if rate else self.BkgCorrected_rawCounts
    
    ax.hist(x = self.mBinX, bins = self.mBin_edges, weights = weights, **kwargs)    
    return
  
  def mpl_errorbar(self, ax, rate = False, bkgCorrected = False, **kwargs):
    weights = self.mBin_rate if rate else self.mBin_rawCounts
    yerr = self.mBin_rateError if rate else np.sqrt(self.mBin_rawCounts)
    if bkgCorrected == True and hasattr(self, "BkgCorrected_rawCounts"):
      weights = self.BkgCorrected_rate if rate else self.BkgCorrected_rawCounts
      yerr = self.BkgCorrected_rateError if rate else self.BkgCorrected_errors
    ax.errorbar(self.mBinX, weights, yerr=yerr, **kwargs)
    return
  
  def GetNP_integral(self, xLow, xHigh, option="rate"):
    """Get integral of numpy histogram between xLow and xHigh."""
    bin_mask = (self.mBinX >= xLow) & (self.mBinX <= xHigh)
    integral = np.sum(self.mBin_rate[bin_mask]) if option == "rate" else np.sum(self.mBin_rawCounts[bin_mask])
    return integral

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