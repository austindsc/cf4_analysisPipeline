## this is now poorly written to only work with LiteRun (from picklelite)

import numpy as np

# this is all the files in a single folder
class RealDataSet(object):

  def __init__(self):
    
    ## only Cs data right now
    self.mCs = []
    self.mCo = []
    
    self.mCoLP = []
    self.mCsLP = []
    
    self.SkipHeights = []
    return
  
  def setup_forComparison(self, Nbins, lowRange, maxRange):
    """Rebin all loaded LiteRun spectra to match MC binning."""
    new_edges = np.linspace(lowRange, maxRange, Nbins + 1)
    for rc in self.mAllRuns:
      if rc.bin_edges is not None and len(rc.bin_centers) != Nbins:
        rc.rebin(new_edges)

  @property
  def mAllRuns(self):
    return self.mCs + self.mCo + self.mCsLP + self.mCoLP

  @property
  def _dataMap(self):
    return {
        'csb': self.mCs,
        'cob': self.mCo,
        'cslp': self.mCsLP,
        'colp': self.mCoLP
    }
  
  @property
  def ValidKeys(self):
    return list(self._dataMap.keys())

  def findRunByHeight(self, height, key="csb", tol=0.1):
    """
    Find a run object by height for a given key.
    
    Args:
        height: The height value to search for
        key: Dataset key ('csb', 'cob', 'cslp', 'colp')
        tol: Tolerance for height matching (default 0.1)
    
    Returns:
        The matching LiteRun/RunClass object, or None if not found
    """
    if key not in self.ValidKeys:
      print(f"ERROR: Invalid key '{key}'. Valid keys are {self.ValidKeys}")
      return None
    
    dataList = self._dataMap.get(key, [])
    if len(dataList) == 0:
      print(f"WARNING: No data loaded for key '{key}'")
      return None
    
    found = next((rc for rc in dataList if np.isclose(rc.mHeight, height, atol=tol)), None)
    if found is None:
      print(f"WARNING: No run found in '{key}' with height {height} (tol={tol})")
    return found

  def print(self, verbose = 0):
    
    print("\n---- RealDataSet print() ----")
    if verbose >= 1:
      if len(self.mCs) > 0:
        print("mCs[0] print():", end = "")
        self.mCs[0].print()
      if len(self.mCo) > 0:
        print("mCo[0] print():", end = "")
        self.mCo[0].print()
      if len(self.mCsLP) > 0:
        print("mCsLP[0] print():", end = "")
        self.mCsLP[0].print()
      if len(self.mCoLP) > 0:
        print("mCoLP[0] print():", end = "")
        self.mCoLP[0].print()
    
    # sketchily print out a table
    wid = 6
    print("\n    h         Nraw        Nsub     rate (Hz)    subRate (Hz)", end = "")
    if len(self.mCs) > 0:
      print("Cs:")
      [print(f"| {rc.mHeight: > {wid}.{1}f} |  {rc.mNEvents:>{wid}.{1}e}  |  {np.sum(rc.BkgCorrected_rawCounts) :>{wid}.{2}e}  |  {rc.GetRate():>{wid}.{2}f} | {rc.GetCorrectedRate():>{wid}.{2}f}") for rc in self.mCs]
    if len(self.mCo) > 0:
      print("\nCo:")
      [print(f"| {rc.mHeight: > {wid}.{1}f} |  {rc.mNEvents:>{wid}.{1}e}  |  {np.sum(rc.BkgCorrected_rawCounts):>{wid}.{2}e}  |  {rc.GetRate():>{wid}.{2}f} | {rc.GetCorrectedRate():>{wid}.{2}f}") for rc in self.mCo]
    if len(self.mCsLP) > 0:
      print("\nCs LP:")
      [print(f"| {rc.mHeight: > {wid}.{1}f} |  {rc.mNEvents:>{wid}.{1}e}  |  {np.sum(rc.BkgCorrected_rawCounts) :>{wid}.{2}e}  |  {rc.GetRate():>{wid}.{2}f} | {rc.GetCorrectedRate():>{wid}.{2}f}") for rc in self.mCsLP]
    if len(self.mCoLP) > 0:
      print("\nCo LP:")
      [print(f"| {rc.mHeight: > {wid}.{1}f} |  {rc.mNEvents:>{wid}.{1}e}  |  {np.sum(rc.BkgCorrected_rawCounts):>{wid}.{2}e}  |  {rc.GetRate():>{wid}.{2}f} | {rc.GetCorrectedRate():>{wid}.{2}f}") for rc in self.mCoLP]
    return
  
  def ModifyHeights(self, hScale, abs = False):
    for rc in self.mAllRuns:
      rc.mHeight = hScale - rc.mHeight
      if abs:
        rc.mHeight = np.abs(rc.mHeight)
    return
  
  # def findRunByHeight(self, height, key = "csb"):
    # # [print(rc.mHeight) for rc in self.gammaMap.get(key)]
    # found = next((rc for rc in self.gammaMap.get(key) if rc.mHeight == height), None)
    # if found == None:
    #   print(f"ERROR: no mc object found in rc {key} with height {height}")
    # return found
  
  
  # def plot_singleRows(self, fig: plt.figure, selection = "mCs",):
  #   Nrows = len(getattr(self, selection))

  #   for i, rc in enumerate(getattr(self, selection)):
  #     ax = fig.add_subplot.hist(Nrows, 1, i+1)
  #     ax.hist(x = rc.mBinX, bins = rc.mBin_edges, weights = rc.mBin_rawCounts, histtype = "step", label = "data " + str(hScale - rcBare_55cm.mHeight) + " cm")
  #   return
