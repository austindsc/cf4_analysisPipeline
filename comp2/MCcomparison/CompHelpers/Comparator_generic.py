"""
Comparator_generic.py - MC to Real Data Comparison Framework

Standalone Functions:
  calc_pvalue(chi2, nBins, nFitParams)  - Calculate p-value from chi2 and degrees of freedom
  calc_reduced_chi2(chi2, nBins, nFitParams) - Calculate reduced chi2

Classes:
  rcMCpair         - Pairs a real data run with an MC run for comparison at a given height
  StatsManager     - Stores chi2, nBins, scaling values
  Comparator       - Main class for managing multiple rcMCpair objects and performing fits

rcMCpair Methods:
  __init__(rc, mc)                   - Initialize pair with real data and MC run objects
  nBins (property)                   - Number of bins used in comparison (from stats)
  BinningInfo (property)             - Dict with binning info + metadata (height, source, chi2)
  BinFromXvalue(xValue)              - Convert x-axis value to bin index
  XFromBinValue(binIndex)            - Convert bin index to x-axis value
  SetBinRanges(firstBin, lastBin, mode, verbose)         - Set analysis bin range (by bin or x-value)
  SetBinRanges_byCountLimit(firstBin, countLimit, option, verbose) - Auto-set lastBin by count threshold
  GetPrefac()                        - Get prefactor (runtime / Nbeam) for scaling MC to data
  GuessScaling()                     - Estimate scaling factor from area ratio
  GetMaxBin(scaling, rawCounts)      - Get max bin value for plotting
  computeChi2(scaling, chiOption)    - Calculate chi2 for this pair; returns (chi2, nBins)

StatsManager Methods:
  __init__()                         - Initialize with default values
  toDict (property)                  - Return stats as a dictionary

Comparator Methods:
  __init__(realData)                 - Initialize with real data set
  UsedPairs (property)               - Flat list of all rcMCpair objects for keysUsed
  UsedPairLists (property)           - List of pair lists, one per key in keysUsed
  ValidKeys (property)               - List of valid key names (csb, cob, cslp, colp)
  loadNewG4(mc, verbose)             - Load a new MC dataset
  ResetPairs()                       - Clear all pair lists
  auto_map_by_common_heights(keys, verbose, heightSkipList, hBuffer) - Auto-create pairs by matching heights
  GetNBins(key)                      - Get total nBins for a single key
  GetTotalNBins(keys)                - Get total nBins across keys (defaults to keysUsed)
  print(mcVerb, dataVerb)            - Print summary of comparator state
  SetBinRanges(firstBin, lastBin, verbose, mode)         - Set bin ranges for all used pairs
  SetBinRanges_byCountLimit(firstBin, countLimit, option, verbose) - Set bin ranges by count limit for all pairs
  SetBinRange_forRun(key, height, firstBin, lastBin, mode, verbose, hBuffer) - Override bin range for one specific run
  GuessOverallScaling()              - Compute initial scaling guess for each key
  GetScaling(key)                    - Get stored scaling for a key
  GetStats(key)                      - Get StatsManager object for a key
  GetPvalue(key, nFitParams)         - Get p-value for a key's chi2 fit
  SetScaling(key, scaling)           - Set scaling for a key
  CalcKeyChi2(key, scaling, chiOption)   - Calculate and store chi2 for a single key
  TotalChi2_singleScaling(scaling, keys, keyRatios, chiOption) - Calc total chi2 with one scaling and ratios
  GetScaledGuess(keys, keyRatios)    - Get scaled initial guesses for dual-source fitting
  GetTotalChi2(keys, chiOption)      - Get total chi2 across keys (defaults to keysUsed)
  minimizeRateScaling(key, statsOption, kwargs, manual, verbose) - Minimize chi2 to find best scaling
"""

import os
import sys
import scipy.optimize as spOpt

import numpy as np
from scipy.stats import chi2 as chi2_dist


########################################
# Standalone statistical functions
########################################

def calc_pvalue(chi2, nBins, nFitParams=1):
  """Calculate p-value for chi2 test. Works with scalars or arrays.
  
  Args:
      chi2: Chi-squared value (scalar or array)
      nBins: Number of bins in comparison (scalar or array)
      nFitParams: Number of fitted parameters (default=1 for scaling).
                  NDOF = nBins - nFitParams
  Returns:
      p-value (probability of observing chi2 this extreme or more)
  """
  chi2 = np.asarray(chi2)
  nBins = np.asarray(nBins)
  ndof = nBins - nFitParams
  
  # Handle invalid values (ndof <= 0 or chi2 < 0)
  result = np.where((ndof > 0) & (chi2 >= 0), chi2_dist.sf(chi2, ndof), np.nan)
  
  # Return scalar if input was scalar
  return result.item() if result.ndim == 0 else result


def calc_reduced_chi2(chi2, nBins, nFitParams=1):
  """Calculate reduced chi2. Works with scalars or arrays.
  
  Args:
      chi2: Chi-squared value (scalar or array)
      nBins: Number of bins in comparison (scalar or array)
      nFitParams: Number of fitted parameters (default=1 for scaling).
                  NDOF = nBins - nFitParams
  Returns:
      Reduced chi2 (chi2 / NDOF)
  """
  chi2 = np.asarray(chi2)
  nBins = np.asarray(nBins)
  ndof = nBins - nFitParams
  
  # Handle invalid values (ndof <= 0)
  result = np.where(ndof > 0, chi2 / ndof, np.nan)
  
  # Return scalar if input was scalar
  return result.item() if result.ndim == 0 else result

def solidAngleFraction(thetaDeg):
  thetaRad = np.deg2rad(thetaDeg)
  return (1 - np.cos(thetaRad)) / 2


########################################
# should have this working as a two way object (rc and mc also know about it)
class rcMCpair(object):

  def __init__(self, rc, mc):
    self.mc = mc
    self.rc = rc
    self.height = mc.mRunPars.gammaHeight
    self.countLimit = 10
    self.Nbeam = mc.Nbeam
    self.stats = StatsManager()  # Each pair has its own stats
    return
  
  @property
  def nBins(self):
    """Convenience property for number of bins in comparison."""
    return self.stats.nBins
  
  @property
  def BinningInfo(self):
    """Return binning info as a dict (includes stats + metadata)."""
    return {
        'firstBin': self.stats.firstBin,
        'lastBin': self.stats.lastBin,
        'nBins': self.stats.nBins,
        'firstBinX': self.rc.mBinX[self.stats.firstBin] if self.stats.nBins > 0 else None,
        'lastBinX': self.rc.mBinX[self.stats.lastBin] if self.stats.nBins > 0 else None,
        'chi2': self.stats.chi2,
        'height': self.height,
        'source': self.mc.mRunPars.gammaSrc,
    }
  
  def BinFromXvalue(self, xValue):
    index = np.abs(self.rc.mBinX - xValue).argmin()
    return index
  
  def XFromBinValue(self, binIndex):
    return self.rc.mBinX[binIndex]
  
  def SetBinRanges(self, firstBin = 1, lastBin = None, mode = "bin", verbose = 0):
    """
    Set analysis bin ranges for data-MC comparison.
    
    Args:
        firstBin: Starting bin index (0-based) or x-axis value, depending on mode
        lastBin: Ending bin index (0-based, inclusive) or x-axis value.
                 If None, the upper edge is chosen by a count-limit criterion.
        mode: "bin" for bin indices, "x" or "value" for x-axis values
        verbose: Verbosity level
        
    Returns:
        dict with keys: 'status', 'firstBin', 'firstBinX', 'lastBin', 'lastBinX', 
                        'nBins', 'countsInLowestBin', 'height', 'source'
    """
    # Use consistent bin definitions already stored on rc and mc:
    #   real data: rc.mBin_Edges, rc.BkgCorrected_rawCounts, rc.mBinX
    #   MC:        mc.mBin_Edges, mc.mRaw_Counts,  mc.mBinX
    if mode == 'bin':
      self.firstBin = firstBin
      self.lastBin = lastBin
    elif mode in ['x', 'value']:
      self.firstBin = self.BinFromXvalue(firstBin)
      self.lastBin = self.BinFromXvalue(lastBin)
    
    # Update stats
    self.stats.firstBin = self.firstBin
    self.stats.lastBin = self.lastBin
    self.stats.nBins = self.lastBin - self.firstBin + 1
    
    # Get counts in the lowest (last) bin for MC
    countsInLowestBin = int(self.mc.mRaw_Counts[self.lastBin])
    
    # Build base result (includes stats info + extra metadata)
    result = {
        'status': 0,
        'firstBin': self.firstBin,
        'firstBinX': self.rc.mBinX[self.firstBin],
        'lastBin': self.lastBin,
        'lastBinX': self.rc.mBinX[self.lastBin],
        'nBins': self.stats.nBins,
        'countsInLowestBin': countsInLowestBin,
        'height': self.height,
        'source': self.mc.mRunPars.gammaSrc,
    }
    
    ## check if any bins within range have low stats in MC
    ## use < (not <=) for consistency with SetBinRanges_byCountLimit which uses >
    belowStat_mask = self.mc.mRaw_Counts[self.firstBin:self.lastBin+1] < self.countLimit
    if np.any(belowStat_mask):
      lowStatBins = np.where(belowStat_mask)[0] + self.firstBin
      lowStatXvalues = self.mc.mBinX[lowStatBins]
      lowStatCounts = self.mc.mRaw_Counts[lowStatBins]
      
      result['status'] = -1
      result['lowStatBins'] = lowStatBins
      result['lowStatXvals'] = lowStatXvalues
      result['lowStatCounts'] = lowStatCounts
      
      if verbose > -1:
        print(f"WARNING in rcmcPair.SetBinRanges: MC counts below {self.countLimit} in bins {lowStatBins}")
        print(f"  x values: {lowStatXvalues}")
        print(f"  MC counts: {lowStatCounts}")
        print(f"  height: {self.height} cm, source: {self.mc.mRunPars.gammaSrc}\n")
      
    return result
  
  def SetBinRanges_byCountLimit(self, firstBin = 1, countLimit = 10, option = "mc", verbose = 1):
    """Set lastBin automatically based on count limit.
    
    Finds the last contiguous bin starting from firstBin where all bins have 
    counts >= countLimit. This ensures no bins within the range fall below the limit.
    
    Args:
        firstBin: Starting bin index (0-based)
        countLimit: Minimum count threshold for bins
        option: Which histogram(s) to check:
            - "mc" or "MC": Use MC counts only
            - "data": Use data counts only  
            - "both": Require BOTH data AND MC to have counts >= countLimit
                      (uses element-wise minimum)
        verbose: Verbosity level
        
    Returns:
        dict with keys: 'status', 'firstBin', 'firstBinX', 'lastBin', 'lastBinX', 
                        'nBins', 'countsInLowestBin', 'height', 'source'
    """
    self.countLimit = countLimit
    self.firstBin = firstBin # by bin, not x value
    
    data_counts = self.rc.BkgCorrected_rawCounts
    if option == "mc" or option == "MC":
      data_counts = self.mc.mRaw_Counts
    elif option == "both":
      # Use element-wise minimum of data and MC counts
      data_counts = np.minimum(self.rc.BkgCorrected_rawCounts, self.mc.mRaw_Counts)
    
    # Find last contiguous bin from firstBin where all bins have counts >= countLimit
    # Start from firstBin and go until we hit a bin below the limit
    belowLimit_mask = data_counts[firstBin:] < countLimit
    has_bad_bin = np.any(belowLimit_mask)
    if has_bad_bin:
      # First bin (relative to firstBin) that falls below limit
      first_bad_bin_relative = np.argmax(belowLimit_mask)
      lastBin = firstBin + first_bad_bin_relative - 1
    else:
      first_bad_bin_relative = None
      # All bins from firstBin onward are above limit
      lastBin = len(data_counts) - 1
    
    # Ensure lastBin is at least firstBin (edge case)
    if lastBin < firstBin:
      lastBin = firstBin
    
    # Call SetBinRanges to set the bins (suppress its output, we handle messaging here)
    binResult = self.SetBinRanges(firstBin = firstBin, lastBin = lastBin, mode = "bin", verbose = -1)
    
    if self.stats.nBins <= 3:
      binResult['status'] = -1
      if verbose >= 0:
        print(f"WARNING in rcmcPair.SetBinRanges_byCountLimit: Only {self.stats.nBins} bins above count limit {countLimit} starting from bin {firstBin}.")
    
    if verbose >= 1:
      src = self.mc.mRunPars.gammaSrc
      h = self.height
      last_x = binResult['lastBinX']
      if not has_bad_bin:
        limit_str = "end-of-array"
        cutoff_pe_str = f"last PE={last_x:.2f}"
        is_mc_limited = False
      else:
        cutoff_abs = firstBin + first_bad_bin_relative
        cutoff_x = self.rc.mBinX[cutoff_abs]
        if option in ("mc", "MC"):
          limit_str = "MC-limited"
          is_mc_limited = True
        elif option == "data":
          limit_str = "data-limited"
          is_mc_limited = False
        else:  # "both" — check which side was actually limiting
          d_low = self.rc.BkgCorrected_rawCounts[cutoff_abs] < countLimit
          m_low = self.mc.mRaw_Counts[cutoff_abs] < countLimit
          if d_low and m_low:
            limit_str = "data+MC-limited"
            is_mc_limited = True
          elif d_low:
            limit_str = "data-limited"
            is_mc_limited = False
          else:
            limit_str = "MC-limited"
            is_mc_limited = True
        cutoff_pe_str = f"cutoff PE={cutoff_x:.2f} (last included PE={last_x:.2f})"
      if verbose == 1 or (verbose >= 2 and is_mc_limited):
        print(f"  h={h:.1f} {src}: {limit_str}, {cutoff_pe_str}, nBins={self.stats.nBins}")
    
    return binResult

  ##----------##
  def GetPrefac(self):
    return self.rc.GetRunTime()/self.Nbeam
    # return self.rc.GetRunTime() / 1.0
  
  def SolidAngleRgamma(self, rgamma):
    angle = float(self.mc.mRunPars.gammaAngle)
    solidFrac = solidAngleFraction(thetaDeg=angle)
    return solidFrac * rgamma
  
  def FixedScaling(self, rgamma):
    return self.SolidAngleRgamma(rgamma) * self.GetPrefac()
  
  def GuessScaling(self): # this is for raw counts in data
    # analytic scaling
    dataArea = sum(self.rc.BkgCorrected_rawCounts[self.firstBin:self.lastBin+1]) # counts
    # xx
    mcAreaByRuntime = np.sum(self.mc.mRaw_Counts[self.firstBin:self.lastBin+1]) * self.GetPrefac() # counts x time / N

    self.mcScale = dataArea / mcAreaByRuntime
    ## all this is confusing
    ## Mdata = M_mc = M_mc_saved * scale * tdata / nbeam
    ## scale = Mdata * nbeam / (M_mc_saved * tdata)
    return self.mcScale # scale is proportional to 1/runtime... kinda. has units 1/time to turn MC to rate

  # def GetMCSpectrum(scaling = scaling, bin0 = 0, ):
  #   return self.mc.mRaw_Counts * scaling * self.rc.GetRunTime() / self.mc.Nbeam

  def GetMaxBin(self, scaling, rawCounts = True):
    if rawCounts:
      a = max(np.maximum(self.rc.BkgCorrected_rawCounts, scaling * self.mc.mRaw_Counts * self.rc.GetRunTime() / self.Nbeam))
      return a
    else:
      a = max(np.maximum(self.rc.BkgCorrected_rawCounts / self.rc.GetRunTime(), scaling * self.mc.mRaw_Counts ))
      return a

  ##----------##
  ## sets the pair's stats values to the chi2 for this pair, and the scaling used
  def computeChi2(self, scaling, chiOption = None): #, rawCounts = True
    if not hasattr(self, "lastBin") or not hasattr(self, "firstBin"):
      print("WARNING in rcmcPair.computeChi2: first and lastBin not set. Call SetBinRanges() first")
      return -1, 0

    if scaling is None or scaling < 0:
      print("something wrong with scaling value in pair.computeChi2")
      return -1, 0
  
    mcBins = scaling * self.GetPrefac() * (self.mc.mRaw_Counts[self.firstBin:self.lastBin+1])
    dataBins = self.rc.BkgCorrected_rawCounts[self.firstBin:self.lastBin+1]
    if chiOption == None:
      chi2_val = np.sum((dataBins - mcBins)**2/mcBins)
      
    ## add other chi2 options here (e.g. with data uncertainty, or using poisson likelihood instead of chi2)
    
    # Update stats
    self.stats.chi2 = chi2_val
    self.stats.nuisScaling = scaling
    
    return self.stats.chi2, self.stats.nBins

########################################
class StatsManager(object):
  def __init__(self):
    self.nuisScaling = None
    self.nuisGuess = 1.0
    self.chi2 = -1
    self.nBins = 0
    self.firstBin = 0
    self.lastBin = 0
    return
  
  def reset(self):
    """Reset all stats to defaults (for reuse across loop iterations)."""
    self.nuisScaling = None
    self.nuisGuess = 1.0
    self.chi2 = -1
    self.nBins = 0
    self.firstBin = 0
    self.lastBin = 0

  @property
  def toDict(self):
    """Return stats as a dictionary for easy access."""
    return {
        'chi2': self.chi2,
        'nBins': self.nBins,
        'firstBin': self.firstBin,
        'lastBin': self.lastBin,
        'nuisScaling': self.nuisScaling,
        'nuisGuess': self.nuisGuess,
    }

########################################
# we use keysUsed to track with datasets are being used (csb, cob, cslp, colp)
# UsedPairs autoupdates depending on what keysUsed is, keysUsed gets set in auto_map_by_common_heights
class Comparator(object):

  def __init__(self, realData):

    self.mMCDataSet = None
    self.mRealData = realData  # this never changes

    # Cesium and Cobalt comparison pairs (generic, not Bare/LP specific)
    self.csPairs = []
    self.csLPpairs = []
    
    self.coPairs = []
    self.coLPpairs = []
    
    # Key-based lookup for pair lists (initialized once, values are refs to actual lists)
    self._pairListMap = {
        'csb': self.csPairs,
        'cob': self.coPairs,
        'cslp': self.csLPpairs,
        'colp': self.coLPpairs
    }
    
    self.keysUsed = [] ## very important?

    self.statsCs = StatsManager()
    self.statsCo = StatsManager()
    self.statsCsLP = StatsManager()
    self.statsCoLP = StatsManager()
    self.statsTotal = StatsManager()

    self.totalRedChi2 = -1
    return
  
  @property
  def UsedPairs(self): # list of RCMC pairs being used
    result = []
    for key in self.keysUsed:
      pairList = self._pairListMap.get(key, [])
      if len(pairList) > 0:
        result.extend(pairList)
    return result
  
  @property
  def UsedPairLists(self): # list of lists of RCMC pairs being used (separate for csb, cob, cslp, colp)
    result = []
    for key in self.keysUsed:
      result.append(self._pairListMap.get(key, []))
    return result
  
  @property
  def _dataMap(self):
    return {
        'csb': self.mRealData.mCs,
        'cob': self.mRealData.mCo,
        'cslp': self.mRealData.mCsLP,
        'colp': self.mRealData.mCoLP
    }
    
  @property
  def _statsMap(self):
    return {
        'csb': self.statsCs,
        'cob': self.statsCo,
        'cslp': self.statsCsLP,
        'colp': self.statsCoLP
    }
  
  @property
  def ValidKeys(self):
    return list(self._pairListMap.keys())
  
  def _validateKeys(self, keys, checkUsed=True):
    """Normalize and validate keys.
    
    Args:
        keys: A single key string, list of keys, or None (defaults to keysUsed).
        checkUsed: If True, also check that keys are a subset of keysUsed.
    Returns:
        List of validated keys.
    Raises:
        ValueError if any key is invalid or (if checkUsed) not in keysUsed.
    """
    if keys is None:
      return list(self.keysUsed)
    if isinstance(keys, str):
      keys = [keys]
    for key in keys:
      if key not in self.ValidKeys:
        raise ValueError(f"Comparator: Invalid key '{key}'. Valid keys are {self.ValidKeys}.")
    if checkUsed and not set(keys).issubset(self.keysUsed):
      missing = set(keys) - set(self.keysUsed)
      raise ValueError(f"Comparator: Keys {missing} not in keysUsed {self.keysUsed}. Call auto_map_by_common_heights() first.")
    return list(keys)

  ### bookkeeping ###
  def loadNewG4(self, mc, verbose = 1):
    # load it
    self.mMCDataSet = mc

    # give it a physPars copy
    return True

  def ResetPairs(self):
    self.csPairs.clear()
    self.coPairs.clear()
    self.csLPpairs.clear()
    self.coLPpairs.clear()
    return

  ## to restore just rerun auto_map_by_common_heights with the same keys and heightSkipList (if any)
  def RemovePairByHeight(self, key, height, hBuffer=0.1):
    """Remove pairs matching a given key and height from the active pair list.

    Args:
        key: Dataset key ('csb', 'cob', 'cslp', 'colp')
        height: Height value to match for removal
        hBuffer: Tolerance for height matching (default 0.1)
    Returns:
        Number of pairs removed.
    """
    if key not in self.ValidKeys:
      raise ValueError(f"Comparator.RemovePairByHeight: Invalid key '{key}'. Valid keys are {self.ValidKeys}.")
    pairList = self._pairListMap[key]
    before = len(pairList)
    pairList[:] = [p for p in pairList if abs(p.height - height) >= hBuffer]
    removed = before - len(pairList)
    print(f"RemovePairByHeight: removed {removed} pair(s) from '{key}' at height ~{height} cm.")
    return removed

  def auto_map_by_common_heights(self, keys=[], verbose = False, heightSkipList=[], hBuffer = 0.1):
    """
    Automatically build rcMCpair mappings by matching heights in MC and real data.
    
    Args:
        cs: If True, create Cs pairs
        co: If True, create Co pairs
        verbose: Print matching details
    """
    for key in keys:
      if key not in self.ValidKeys:
        raise ValueError(f"invalid key passed")
    self.keysUsed = keys

    if self.mMCDataSet is None:
      raise RuntimeError("Comparator: MC dataset not loaded. Call loadNewG4() first.")
    
    if self.mRealData is None:
      raise RuntimeError("Comparator: Real data not loaded.")

    ## doing everything below generically
    ## _pairListMap is a dict initialized in __init__ whose values are refs to the actual lists
    dataMap = self._dataMap
    
    for givenKey in keys:
      pairList = self._pairListMap[givenKey]  # Get the actual list reference
      pairList.clear()  # Reset the actual list
      for rc in dataMap[givenKey]:
        if not hasattr(rc, "mHeight"):
          raise RuntimeError(f"RealDataSet object missing mHeight attribute for {givenKey} data.")
        if np.any(np.isclose(rc.mHeight, heightSkipList, atol=hBuffer)):
          if verbose:
            print(f"Skipping {givenKey} data at height {rc.mHeight} cm as it is in the skip list.")
          continue
        matching_mc = self.mMCDataSet.findRunByHeight(rc.mHeight, key=givenKey)
        if matching_mc is not None:
          pairList.append(rcMCpair(rc, matching_mc))  # Append to the actual list
          if verbose:
            print(f"Mapped {givenKey} data at height {rc.mHeight} cm to MC at height {matching_mc.mRunPars.gammaHeight} cm.")
        else:
          if verbose:
            print(f"No matching MC found for {givenKey} data at height {rc.mHeight} cm.")
    return True
  
  def GetNBins(self, key = "csb"):
    """Get total nBins for a single key."""
    pairList = self._pairListMap.get(key, [])
    return sum(pair.stats.nBins for pair in pairList)
  
  def GetTotalNBins(self, keys=None):
    """Get total nBins across specified keys.
    
    Args:
        keys: List of keys to sum over. If None, uses keysUsed.
    Returns:
        Total number of bins across all pairs for the given keys.
    """
    if keys is None:
      keys = self.keysUsed
    total = 0
    for key in keys:
      total += self.GetNBins(key)
    return total
  
  ### Functionality prints ###
  def print(self, mcVerb = True, dataVerb = True): # incomplete
    print("\n---- COMPARATOR.print() ----")
    if mcVerb or dataVerb:
      print("\t--- subcalled ---")
      self.mMCDataSet.print() if mcVerb else None
      self.mRealData.print(verbose = dataVerb, hScale = self.hScale) if dataVerb else None

    if len(self.keysUsed) == 0:
      print("No keys used yet in comparator, so no pairs to print.\n")
      return
    for key in self.keysUsed:
      pairList = self._pairListMap.get(key, [])
      print(f"Heights for {key} w {len(pairList)} entries:")
      [print("{", self.hScale - pair.rc.mHeight, " : ", pair.mc.mRunPars.gammaHeight, "}", sep = "", end = "\n") for pair in pairList]
    return
    
  ### Statistics and fitting ###
  def SetBinRanges(self, firstBin = 1, lastBin = 10, verbose = 1, mode = "value"):
    """
    Set bin ranges for all used pairs.
    
    Returns:
        list of result dicts from each rcMCpair.SetBinRanges() call.
        Each dict has: 'status', 'firstBin', 'lastBin', 'firstX', 'lastX', 'nBins', 
                       'height', 'source', 'run', and if status=-1: 'lowStatBins', 
                       'lowStatXvals', 'lowStatCounts'
    """
    results = []
    for rcmc in self.UsedPairs:
      res = rcmc.SetBinRanges(firstBin=firstBin, lastBin=lastBin, verbose=verbose, mode=mode)
      results.append(res)
    
    nIssues = sum(1 for r in results if r['status'] != 0)
    if nIssues > 0 and verbose > -1:
      print(f"WARNING in Comparator.SetBinRanges: {nIssues} pairs had issues with bin range setting.")
    elif verbose > -1:
      print(f"---\nComparator.SetBinRanges: all {len(results)} pairs set successfully within range {firstBin} to {lastBin} ({mode} mode).\n---")
    return results
  
  def SetBinRanges_byCountLimit(self, firstBin = 1, countLimit = 10, option = "mc", verbose = 1):
    """Set bin ranges for all used pairs based on count limit.
    
    Returns:
        list of result dicts from each rcMCpair.SetBinRanges_byCountLimit() call.
    """
    if verbose >= 1:
      print(f"Comparator.SetBinRanges_byCountLimit: {len(self.UsedPairs)} pairs, countLimit={countLimit}, option='{option}'")
    results = []
    for rcmc in self.UsedPairs:
      res = rcmc.SetBinRanges_byCountLimit(firstBin = firstBin, countLimit = countLimit, option = option, verbose = verbose)
      results.append(res)
    if verbose >= 1:
      nIssues = sum(1 for r in results if r['status'] != 0)
      print(f"  -> {len(results)} pairs set; {nIssues} with status != 0.")
    return results
  
  def SetBinRange_forRun(self, key, height, firstBin=None, lastBin=None, mode="bin", verbose=1, hBuffer=0.1):
    """Override bin range for a single pair identified by key and height.

    Call this after SetBinRanges_byCountLimit to patch one specific run.

    Args:
        key:      Key string (e.g. "cslp", "colp").
        height:   Height of the run to target (cm).
        firstBin: New lower bound (bin index or x-value depending on mode).
                  If provided without lastBin, only the lower boundary is changed;
                  the existing lastBin (already set by the global SetBinRanges call)
                  is preserved unchanged.
        lastBin:  New upper bound. If provided together with firstBin, sets both
                  directly.
        mode:     "bin" or "value"/"x" — passed to SetBinRanges when lastBin given.
        verbose:  Verbosity level.
        hBuffer:  Tolerance for height matching (cm).

    Returns:
        Result dict from rcMCpair.SetBinRanges or SetBinRanges_byCountLimit.
    Raises:
        ValueError: if no pair is found for the given key/height.
    """
    if key not in self.ValidKeys:
      raise ValueError(f"SetBinRange_forRun: invalid key '{key}'.")
    pairList = self._pairListMap.get(key, [])
    target = None
    for pair in pairList:
      if np.isclose(pair.height, height, atol=hBuffer):
        target = pair
        break
    if target is None:
      raise ValueError(
        f"SetBinRange_forRun: no pair found for key='{key}' at height={height} cm "
        f"(buffer={hBuffer}). Available heights: {[p.height for p in pairList]}"
      )

    if firstBin is not None and lastBin is not None:
      result = target.SetBinRanges(firstBin=firstBin, lastBin=lastBin, mode=mode, verbose=verbose)
    elif firstBin is not None:
      # Only override the lower boundary; preserve the upper boundary already set on this pair.
      # target.lastBin is always stored as a bin index, so normalize firstBin to a bin
      # index as well (handles mode="value") and then call SetBinRanges in "bin" mode.
      if not hasattr(target, 'lastBin'):
        raise RuntimeError(
          f"SetBinRange_forRun: pair at key='{key}' h={height}cm has no lastBin yet. "
          "Call SetBinRanges_byCountLimit (or SetBinRanges) globally first."
        )
      if mode in ('x', 'value'):
        firstBinIdx = target.BinFromXvalue(firstBin)
      else:
        firstBinIdx = firstBin
      result = target.SetBinRanges(firstBin=firstBinIdx, lastBin=target.lastBin, mode="bin", verbose=verbose)
    elif lastBin is not None:
      result = target.SetBinRanges(firstBin=target.firstBin, lastBin=lastBin, mode=mode, verbose=verbose)
    else:
      raise ValueError("SetBinRange_forRun: provide at least one of firstBin or lastBin.")

    if verbose >= 1:
      print(f"SetBinRange_forRun: patched key='{key}' h={height}cm -> "
            f"bins [{result['firstBin']}, {result['lastBin']}] "
            f"(x: {result['firstBinX']:.2f} to {result['lastBinX']:.2f}), "
            f"nBins={result['nBins']}")
    return result

  def GuessOverallScaling(self, keys=None): # call guess scaling for each pair to get overall scaling
    """Compute initial scaling guess for specified keys.
    
    Args:
        keys: A single key string, list of keys, or None (defaults to keysUsed).
    """
    keys = self._validateKeys(keys, checkUsed=True)
    
    for key in keys:
      pairList = self._pairListMap.get(key)
      statsObj = self._statsMap.get(key)
      
      statsObj.nuisGuess = np.sum([rcmc.GuessScaling() for rcmc in pairList])/len(pairList)
    return True

  def GetScaling(self, key):
    statsObj = self._statsMap.get(key)
    return statsObj.nuisScaling
  
  def GetStats(self, key):
    statsObj = self._statsMap.get(key)
    return statsObj
  
  def GetPvalue(self, key, nFitParams=1):
    """Get p-value for a specific key's chi2 fit.
    
    Args:
        key: Dataset key ("colp", "cslp", etc.)
        nFitParams: Number of fitted parameters (default=1 for scaling)
    Returns:
        p-value
    """
    stats = self._statsMap.get(key)
    return calc_pvalue(stats.chi2, stats.nBins, nFitParams)
  
  def SetScaling(self, key, scaling):
    statsObj = self._statsMap.get(key)
    statsObj.nuisScaling = scaling
    return True
  
  # calculates chi2 and nBins for a list (dataset) of rcmc pairs
  # one stats object per pairList (dataset)
  # one stats object per total (key, all pairs in that key)
  def CalcKeyChi2(self, key, scaling=None, chiOption=None):
    pairList = self._pairListMap.get(key)
    stats = self._statsMap.get(key) # this is the global stats object
      
    keyPairValues = np.array([pair.computeChi2(scaling, chiOption=chiOption) for pair in pairList])
    # [[chi2, nBins],...]
    
    self.SetScaling(key, scaling)
    stats.chi2 = np.sum(keyPairValues[:, 0])
    stats.nBins = np.sum(keyPairValues[:, 1])

    return stats.chi2

  ## make one of keyratio 1
  def TotalChi2_singleScaling(self, scaling=None, keys = [], keyRatios=[], chiOption=None):
    totalChi2 = 0
    totalNBins = 0

    for key, ratio in zip(keys, keyRatios):
      self.CalcKeyChi2(key=key, scaling=scaling*ratio, chiOption=chiOption)
      totalChi2 += self._statsMap.get(key).chi2
      totalNBins += self._statsMap.get(key).nBins
    
    self.statsTotal.chi2 = totalChi2
    self.statsTotal.nBins = totalNBins
    self.statsTotal.nuisScaling = scaling
    return totalChi2
  
  ## divide the guesses by the ratios.
  def GetScaledGuess(self, keys, keyRatios=[]):
    scaledGuesses = []
    for key, ratio in zip(keys, keyRatios):
      statsObj = self._statsMap.get(key)
      scaledGuesses.append(statsObj.nuisGuess/ratio)
    return scaledGuesses

  def GetTotalChi2(self, keys=None, chiOption=None):
    """Get total chi2 across specified keys.
    
    Args:
        keys: List of keys to sum over. If None, uses keysUsed.
        chiOption: Chi2 calculation option (passed to CalcKeyChi2)
    Returns:
        Total chi2 across all pairs for the given keys.
    """
    if keys is None:
      keys = self.keysUsed
    
    totalChi2 = 0
    totalNBins = 0

    for key in keys:
      statsObj = self._statsMap.get(key)
      self.CalcKeyChi2(key=key, scaling=statsObj.nuisScaling, chiOption=chiOption)
      
      totalChi2 += statsObj.chi2
      totalNBins += statsObj.nBins
    
    self.statsTotal.chi2 = totalChi2
    self.statsTotal.nBins = totalNBins
    return totalChi2
  
  def minimizeRateScaling(self, key = "co", statsOption = "single", kwargs=None, manual = False, verbose = 0): # nuisance parameter for these dataSets
    # options are same as 'keysUsed' entries: "co", "cs", "colp", "cslp"

    # tryAgain: add a couple attempts with improved minimizer, also add option to do manually
    # (find min of some 1d space like above), also allow option to specify which function gets passed.
    minOptions = {'maxiter': 1000, 'disp': False, 'xatol': 1e-6, 'fatol':1e-6 }
    
    self.GuessOverallScaling(keys=[key]) # sets the guess in the stats object for that key, which is used as x0 for the minimizer
    if statsOption == "single":
      guess = self._statsMap.get(key).nuisGuess
      funToUse = lambda s: self.CalcKeyChi2(key=key, scaling=s)
    
    ## key should be array for dual source
    # if statsOption == "dualSource":
    #   krs = kwargs.get('keyRatios', [1]*len(key))
    #   funToUse = lambda s: self.TotalChi2_singleScaling(scaling=s, keys=key, keyRatios=krs)
    #   guess = np.average(self.GetScaledGuess(keys=key, keyRatios=krs))
      
    keyList = self._validateKeys(key, checkUsed=True)
    
    if manual:
      # result = spOpt.brute(fun = funToUse, ranges = [(0.1*guess, 10.*guess)], Ns = 100)
      scaleRange = np.linspace(0.1*guess, 10.*guess, 400)
      chiValues = [funToUse(s) for s in scaleRange]
      argMin = np.argmin(chiValues)
      funToUse(scaleRange[argMin]) # reset the internals that track nuisance params
      return scaleRange[argMin], chiValues[argMin]
    elif not manual:
      result = spOpt.minimize(fun = funToUse, method = "nelder-mead", options=minOptions, x0 = [guess], bounds = [[0.1*guess, 10.*guess]])
      
    # run again to store min values
    if result.success == True:
      funToUse(result.x[0]) # reset the internals that track nuisance params
      if verbose == 1:
        print(result)
      return result.x[0], result
    else:
      print(result)
      print(f"minimize rate scaling failed for xComp1, key {key}, option {statsOption}")
      return 0, result

  def FixedScaling_Keychi2(self, key, rgamma, chiOption=None):
    pairList = self._pairListMap.get(key)
    stats = self._statsMap.get(key) # this is the global stats object
      
    keyPairValues = np.array([pair.computeChi2(pair.SolidAngleRgamma(rgamma), chiOption=chiOption) for pair in pairList])
    # [[chi2, nBins],...]
    
    stats.chi2 = np.sum(keyPairValues[:, 0])
    stats.nBins = np.sum(keyPairValues[:, 1])
    return stats.chi2
  
  def GetFixedScaling_chi2(self, keys, rgammas, chiOption=None):
    totalChi2 = 0
    totalNBins = 0

    for key, rgamma in zip(keys, rgammas):
      totalChi2 += self.FixedScaling_Keychi2(key=key, rgamma=rgamma, chiOption=chiOption)
      statsObj = self._statsMap.get(key)
      totalChi2 += statsObj.chi2
      totalNBins += statsObj.nBins
    
    self.statsTotal.chi2 = totalChi2
    self.statsTotal.nBins = totalNBins
    return totalChi2
  
  
# find the scaling nuisance parameter (minimize chi2) for the Cs dataset.
# only Bare for now (only have BareMC data)
# def scaleCs(startValue):

# comparator.mc1.mCsBareList[0].mPhysPars.ly
