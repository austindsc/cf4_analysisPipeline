# this will take in a folder and create the MC dataset
# # and an set 
import sys
import os

# Geant4 code and helpers
sys.path.append(os.environ.get("CF4G4"))
import Batch.helpers.FileGetters as FG

from Batch.helpers.G4_runClasses import *
from Batch.helpers.G4_macHelpers import *

# this is all the files in a single folder
class MCdataSet(object):

  def __init__(self, folderPath: str, verbose = 1):
    
    files = FG.get_all_root_files(folderPath)
    self.fpath = folderPath
    if verbose == 1:
      print(files)
    # open everything as FullMCRun
    # sort by mRunPars

    # all runs in one array
    self.mAllRuns = [FullMCRun(macName = f[:-len(".root")], verbose = verbose) for f in files]

    # gamma Map and sorting
    self.buildOutDataArrays()
    return
  
  ########## Management ##########
  # this builds or resets the arrays that sort all MC stuff from mAllRuns
  def buildOutDataArrays(self):

    ## sorting Cs data
    self.mCsBareList = [mc for mc in self.mAllRuns if (mc.mRunPars.lpFilter == False and mc.mRunPars.gammaSrc == "Cs137")]
    self.mCsBareList = sorted(self.mCsBareList, key=lambda rc: rc.mRunPars.gammaHeight) # sorting
    self.mCsLpList = [mc for mc in self.mAllRuns if (mc.mRunPars.lpFilter == True and mc.mRunPars.gammaSrc == "Cs137")]
    self.mCsLpList = sorted(self.mCsLpList, key=lambda rc: rc.mRunPars.gammaHeight) # sorting
    # [print(type(mc.mRunPars.gammaHeight)) for mc in self.mCsBareList]

    ## sorting Co data
    self.mCoBareList = [mc for mc in self.mAllRuns if (mc.mRunPars.lpFilter == False and mc.mRunPars.gammaSrc == "Co57")]
    self.mCoBareList = sorted(self.mCoBareList, key = lambda rc: rc.mRunPars.gammaHeight)
    self.mCoLPList = [mc for mc in self.mAllRuns if (mc.mRunPars.lpFilter == True and mc.mRunPars.gammaSrc == "Co57")]
    self.mCoLPList = sorted(self.mCoLPList, key = lambda rc: rc.mRunPars.gammaHeight)

    self.gammaMap = {"csb": self.mCsBareList, "cslp": self.mCsLpList, "cob": self.mCoBareList, "colp": self.mCoLPList}

    ## other data\
    return
  
  # remove a height from the data lists. Not needed with new mapping strategy
  def popHeights(self, csBare = [], coBare = [], csLP = [], coLP = []):
    self.mCsBareList = [mc for mc in self.mCsBareList if mc.mRunPars.gammaHeight not in csBare]
    self.mCoBareList = [mc for mc in self.mCoBareList if mc.mRunPars.gammaHeight not in coBare]
    return
  
  # important for searching in self for given mc run
  def findRunByHeight(self, height, key = "csb"):
    found = next((mc for mc in self.gammaMap.get(key) if mc.mRunPars.gammaHeight == height), None)
    if found == None:
      print(f"no mc object found in mc {key}:{self.gammaMap.get(key)} with height {height}")
    return found

  ########## Corrections and setup ##########
  def doCAcorrection_BinomialMC(self, nCA):
    done = [mc.doCAcorrection_BinomialMC(nCA) for mc in self.mAllRuns]
    return
  
  def doCAcorrection_GeoMC(self, nCA):
    done = [mc.doCAcorrection_GeometricMC(nCA) for mc in self.mAllRuns]
    return
  
  # for "mNPE", "mNPE_corr", call one of the doCAcorrection first
  def setup_forComparison(self, Nbins = 100, lowRange = 0.5, maxRange = 100.5, variable = "mNPE"):
    done = [mc.setup_forComparison(Nbins, lowRange, maxRange, variable) for mc in self.mAllRuns]
    return
  
  # def compare_dicts(dict1, dict2):
  #   return all(dict1.get(key) == dict2.get(key) for key in dict1.keys())

  def print(self):
    print("\n---- mcDataSet.print() ----")
    print(self.fpath, "\n")
    print("self.mCsBareList[0]", self.mCsBareList[0].mPhysPars)
    print("self.mCsBareList[0]", self.mCsBareList[0].mRunPars)
    # print("Cs Bare heights")
    [print(mc.mRunPars.gammaHeight, end = " ") for mc in self.mCsBareList];

    print("\nself.mCoBareList[0]", self.mCoBareList[0].mRunPars)
    [print(mc.mRunPars.gammaHeight, end = " ") for mc in self.mCoBareList]
    print()
    return
  
  # check that all physPars are the same (they should be)
  # detPars aren't implemented yet (only PDE)
  def check_physPars_andDetPars(self, verbose = 1):

    physParList = [mc.mPhysPars for mc in self.mAllRuns]
    if verbose == 2:
      [print(mc.mPhysPars.__dict__, "\n") for mc in self.mAllRuns]
    a = all(physParList[0] == p1 for p1 in physParList[1:])
    if a == True:
      if verbose:
        print("\nphys Pars consistent for", self.fpath)
        print(self.mCsBareList[0].mPhysPars)
      return a
    else:
      print("ERROR, physPars inconsistent for", self.fpath)
    return a