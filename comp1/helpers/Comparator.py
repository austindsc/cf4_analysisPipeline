import os
import sys
import scipy.optimize as spOpt
import json5

import numpy as np

##  tools for real data
sys.path.append(os.environ.get("CF4HOME")+"/code1Analysis/")
sys.path.append(os.environ.get("CF4HOME"))

# from helpers.anaStructs import * # includes SORTING class/namespace only needed in old do_mapping

# should have this working as a two way object (rc and mc also know about it)
class rcMCpair(object):

  def __init__(self, rc, mc):
    self.mc = mc
    self.rc = rc
    self.height = mc.mRunPars.gammaHeight
    self.countLimit = 10
    self.ndof = -1
    self.Nbeam = mc.Nbeam
    return
  
  def SetBinRanges(self, firstBin = 1, verbose = 0):

    self.firstBin = firstBin

    find_index = lambda arr, X: next((i for i, value in enumerate(arr) if value < X), len(arr) - 1)
    dataLastBin = find_index(self.rc.mBin_rawCounts, self.countLimit)
    mcLastBin = find_index(self.mc.mRaw_Counts, self.countLimit)
    self.lastBin = min(dataLastBin, mcLastBin)

    # safety. shouldn't enter this?
    while self.mc.mRaw_Counts[self.lastBin] < self.countLimit or self.rc.mBin_rawCounts[self.lastBin] < self.countLimit:
      self.lastBin -= 1

    self.ndof = self.lastBin - self.firstBin

    if mcLastBin < dataLastBin:
      print("YOO data limited in:", self.mc.mRunPars.fname) if verbose else None
      print("   mc : data lastBin count = ", self.mc.mRaw_Counts[self.lastBin+1], " : ", self.rc.mBin_rawCounts[self.lastBin+1]) if verbose else None
      return self.mc.mRunPars.fname
    return 1

  ##----------##
  def GuessScaling(self): # this is for raw counts
    # analytic scaling
    dataArea = sum(self.rc.mBin_rawCounts[self.firstBin:self.lastBin])
    # xx
    mcAreaByRuntime = np.sum(self.mc.mRaw_Counts[self.firstBin:self.lastBin]) * self.rc.GetRunTime() / self.Nbeam

    self.mcScale = dataArea / mcAreaByRuntime
    return self.mcScale # this is now proportional to the runtime. has units 1/time to turn MC to rate

  # def GetMCSpectrum(scaling = scaling, bin0 = 0, ):
  #   return self.mc.mRaw_Counts * scaling * self.rc.GetRunTime() / self.mc.Nbeam

  def GetMaxBin(self, scaling, rawCounts = True):
    if rawCounts:
      a = max(np.maximum(self.rc.mBin_rawCounts, scaling * self.mc.mRaw_Counts * self.rc.GetRunTime() / self.Nbeam))
      return a
    else:
      a = max(np.maximum(self.rc.mBin_rawCounts / self.rc.GetRunTime(), scaling * self.mc.mRaw_Counts ))
      return a

  ##----------##
  def computeChi2_reduced(self, scaling): #, rawCounts = True
    if not hasattr(self, "lastBin") or not hasattr(self, "firstBin"):
      print("WARNING in rcmcPair.computeChi2_reduced: first and lastBin not set. Call SetBinRanges() first")
      return -1

    if scaling is None or scaling < 0:
      print("something wrong with scaling value in pair.computeChi2_reduced")
      return -1
  
    mcBins = scaling * (self.mc.mRaw_Counts[self.firstBin:self.lastBin]) * self.rc.GetRunTime() / self.Nbeam
    dataBins = self.rc.mBin_rawCounts[self.firstBin:self.lastBin]
    self.chi2 = np.sum((dataBins - mcBins)**2/mcBins)
    self.nbinsUsed = self.lastBin - self.firstBin
    self.redChi2 = self.chi2/self.ndof
    return self.chi2, self.ndof

# do i make another object to track chi2 scores, nuisance pars, etc?

########################################
class Comparator(object):

  def __init__(self, realData):

    self.mMCDataSet = None
    self.mRealData = realData # this never changes
    self.loadConfig()
    self.hScale = -1.0

    ## Cesium
    self.csBareMap = [None] * len(self.mRealData.mCsBare)
    # self.csGuess = None
    self.csBareNuis = None
    self.csBareChi2 = -1
    self.csBareNDOF = 0

    ## Co
    self.coBareMap = [None] * len(self.mRealData.mCoBare)
    # self.coGuess = None
    self.coBareNuis = None
    self.coBareChi2 = -1
    self.coBareNDOF = 0

    self.totalRedChi2 = -1
    return
  
  def loadConfig(self, path = os.environ.get("CF4HOME")+"/comp1/" , fname = "config.json"):
    with open(path+fname, 'r') as file:
      self.config = json5.load(file)
    return self.config

  def loadNewG4(self, mc, verbose = 1):
    # load it
    self.mMCDataSet = mc

    # give it a physPars copy
    return True

  # # how to track heights and map them to one another. 
  # def do_mapping(self):

  #   # check lengths
  #   if (len(self.mMCDataSet.mCsBareList) != len(self.mRealData.mCsBare)):
  #     print("mc csBare list not equal to realData list")
  #     return False

  #   # BareCs mapping
  #   for i, mc in enumerate(self.mMCDataSet.mCsBareList): # here mc is a FullMCRun object
  #     h = mc.mRunPars.gammaHeight # get height
  #     rc_withHeight = SORTING.findRCbyHeight(self.hScale - h, self.mRealData.mCsBare) # get real data
  #     pair = rcMCpair(rc_withHeight, mc) # make pair object to map them together
  #     # print(hScale - rc_withHeight.mHeight, pair.mc.mRunPars.gammaHeight == self.hScale - rc_withHeight.mHeight) 
  #     self.csBareMap[i] = pair # assign pair to the map list
    
  #   # doing BareCo mapping
  #   for i, mc in enumerate(self.mMCDataSet.mCoBareList): # here mc is a FullMCRun object
  #     h = mc.mRunPars.gammaHeight # get height
  #     # print(mc.mRunPars)
  #     # print(h)
  #     rc_withHeight = SORTING.findRCbyHeight(self.hScale - h, self.mRealData.mCoBare) # get real data
  #     # print(h, type(rc_withHeight))
  #     pair = rcMCpair(rc_withHeight, mc) # make pair object to map them together
  #     # print(hScale - rc_withHeight.mHeight, pair.mc.mRunPars.gammaHeight == self.hScale - rc_withHeight.mHeight) 
  #     self.coBareMap[i] = pair # assign pair to the map list
    
  #   # print(len(self.csBareMap))
  #   return True
  
  # should also have option to just pass a 'pure' height list
  def do_mapping_fromJson(self, configName = "config.json", subtract = 0.0, verbose = False):
    print("in do_mapping_fromJson") if verbose else None 

    # reload the config file
    self.loadConfig(fname = configName)

    # csBare mapping
    dic_csBare = self.config.get("global_analysis").get("csBareMap")
    self.csBareMap = [None] * len(dic_csBare)
    print("CsBare") if verbose else None
    for i, (mch, datah) in enumerate(dic_csBare.items()):
      mc = self.mMCDataSet.findRunByHeight(float(mch), key = "csb")
      rc = self.mRealData.findRunByHeight(datah, key = "csb")
      if mc != None and rc != None:
        pair = rcMCpair(rc, mc) # make pair object to map them together
        self.csBareMap[i] = pair # assign pair to the map list
        print(mc.mRunPars.gammaHeight, rc.mHeight) if verbose else None
    
    # coBare mapping
    dic_coBare = self.config.get("global_analysis").get("coBareMap")
    self.coBareMap = [None] * len(dic_coBare)
    print("\nCoBare") if verbose else None
    for i, (mch, datah) in enumerate(dic_coBare.items()):
      mc = self.mMCDataSet.findRunByHeight(float(mch), key = "cob")
      rc = self.mRealData.findRunByHeight(datah, key = "cob")
      if mc != None and rc != None:
        pair = rcMCpair(rc, mc) # make pair object to map them together
        self.coBareMap[i] = pair # assign pair to the map list
        print(mc.mRunPars.gammaHeight, rc.mHeight) if verbose else None
    
    if None in self.coBareMap or None in self.csBareMap:
      print("some None elements exist in COMPARATOR object maps")
      return False
    return True

  def print(self, mcVerb = True, dataVerb = True): # incomplete
    print("---- comparator.print() subcalls ----\n")

    self.mMCDataSet.print() if mcVerb else None

    self.mRealData.print(verbose = dataVerb, hScale = self.hScale)

    print("\n---- COMPARATOR.print() ----")
    print("heights in pairMaps:")
    print("csBare:", end = " ")
    [print("{", self.hScale - pair.rc.mHeight, " : ", pair.mc.mRunPars.gammaHeight, "}", sep = "", end = " ")for pair in self.csBareMap]

    print("\ncoBare:", end = " ")
    [print("{", self.hScale - pair.rc.mHeight, " : ", pair.mc.mRunPars.gammaHeight, "}", sep = "", end = " ")for pair in self.coBareMap]
    print("\n")
    return
    
  def SetBinRanges(self, firstBin = 1, verbose = 0):
    a = [rcmc.SetBinRanges(firstBin = firstBin, verbose = verbose) for rcmc in self.csBareMap] # Cs
    b = [rcmc.SetBinRanges(firstBin = firstBin, verbose = verbose) for rcmc in self.coBareMap] # Co

    finalRes = a + b
    runStrings = [run for run in finalRes if isinstance(run, str)]
    if len(runStrings) > 0:
      print("mc stats limited in", runStrings) if verbose > 1 else None
      return runStrings
    return None
  
  def checkNDOF(self):
    return all(p.ndof != 0 for p in self.csBareMap + self.coBareMap)
  
  def GuessOverallScaling(self): # call guess scaling for each pair to get overall scaling
    ## Cs (just bare or both?)
    self.csGuess = np.sum([pair.GuessScaling() for pair in self.csBareMap])/len(self.csBareMap)
    # self.csBareNuis = self.csGuess # this messes up other succesful minimzations

    ## Co (just bare or both?)
    self.coGuess = np.sum([pair.GuessScaling() for pair in self.coBareMap])/len(self.coBareMap)
    # self.coBareNuis = self.coGuess
    return True

  #### 
  # use for Cs nuisance param
  def getCsBarechi2(self, scaling = None):
    if scaling is None:
      scaling = self.csGuess
    else:
      self.csBareNuis = scaling
    
    csBarePairValues = np.array([pair.computeChi2_reduced(scaling) for pair in self.csBareMap]) # do calculation
    self.csBareChi2 = np.sum(csBarePairValues[:, 0])
    self.csBareNDOF = np.sum(csBarePairValues[:, 1])
    self.csBareRedChi2 = self.csBareChi2/(self.csBareNDOF - 1)
    return self.csBareChi2
  
  # use for Co nuisance param
  def getCoBarechi2(self, scaling = None):
    if scaling is None:
      scaling = self.coGuess
    else:
      self.coBareNuis = scaling
      
    coBarePairValues = np.array([pair.computeChi2_reduced(scaling) for pair in self.coBareMap]) # do calculation
    self.coBareChi2 = np.sum(coBarePairValues[:, 0])
    self.coBareNDOF = np.sum(coBarePairValues[:, 1])
    self.coBareRedChi2 = self.coBareChi2/(self.coBareNDOF - 1)
    # print("full Co X2 = ", self.coBareChi2)
    return self.coBareChi2

  ## this one does everything, but nuisance params don't use this exact output
  def getFullReducedChi2(self, useGuess = False):
    if not hasattr(self, "csGuess") or not hasattr(self, "coGuess"):
      self.GuessOverallScaling()
    if useGuess == True: # this also gets called if self.csBareNuis are None, by the init values of csBareNuis as None
      self.getCsBarechi2()
      self.getCoBarechi2()
    else:
      self.getCsBarechi2(scaling = self.csBareNuis)
      self.getCoBarechi2(scaling = self.coBareNuis)

    self.totalRedChi2 = (self.csBareChi2 + self.coBareChi2)/(self.csBareNDOF + self.coBareNDOF - 2)
    return self.totalRedChi2
  
  def minimizeRateScaling(self, option = "CsBare", manual = False, verbose = 0): # nuisance parameter for these dataSets
    # options are CsBare, CoBare, ...
    opts = ["CoBare", "CsBare"]

    # tryAgain: add a couple attempts with improved minimizer, also add option to do manually
    # (find min of some 1d space like above), also allow option to specify which function gets passed.
    minOptions = {'maxiter': 1000, 'disp': False, 'xatol': 1e-6, 'fatol':1e-6 }
    
    self.GuessOverallScaling()

    # Cesium Bare by default
    funToUse = self.getCsBarechi2
    guess = self.csGuess
    result = None
    if option not in opts:
      print("ERROR, wrong option given in minimizeRateScaling")
      return 0
    elif option == "CoBare":
      funToUse = self.getCoBarechi2
      guess = self.coGuess
    
    result = spOpt.minimize(fun = funToUse, method = "nelder-mead", options=minOptions, x0 = [guess], bounds = [[0.1*guess, 10.*guess]])
    
    # run again to store min values
    if result.success == True:
      funToUse(result.x[0]) # reset the internals that track nuisance params (resets self.coBareNuis to track results)
      if verbose == 1:
        print(result)
      return 1, result
    else:
      print(result)
      print(f"minimize rate scaling failed for xComp1, option {option}")
      return 0, result

# find the scaling nuisance parameter (minimize chi2) for the Cs dataset.
# only Bare for now (only have BareMC data)
# def scaleCs(startValue):


# comparator.mc1.mCsBareList[0].mPhysPars.ly