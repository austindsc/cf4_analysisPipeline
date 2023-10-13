# this is all the files in a single folder
class RealDataSet(object):

  def __init__(self):
    
    ## only Cs data right now
    self.mCsBare = []
    self.mCsLP = []

    self.mCoBare = []
    self.mCoLP = []

    return

  def print(self, hSub = 0, verbose = 0):
    
    print("\n---- RealDataSet print() ----") if verbose > 0 else None
    if verbose > 1:
      print("mCsBare[0] print():", end = "")
      self.mCsBare[0].print()

      print("\nmCoBare[0] print()", end = "")
      self.mCoBare[0].print()
    
    if verbose > 0:
      wid = 6
      print("    h         Nraw        Nsub     rate (Hz)    subRate (Hz)")
      print("csBare:")
      [print(f"| {hSub - rc.mHeight: > {wid}.{1}f} |  {rc.tt.GetEntries():>{wid}.{1}e}  |  {rc.GetCorrSpectrum().GetEntries():>{wid}.{2}e}  |  {rc.GetRate():>{wid}.{2}f} | {rc.mCorrBkgRate:>{wid}.{2}f}") for rc in self.mCsBare]

      print("\ncoBare:")
      [print(f"| {hSub - rc.mHeight: > {wid}.{1}f} |  {rc.tt.GetEntries():>{wid}.{1}e}  |  {rc.GetCorrSpectrum().GetEntries():>{wid}.{2}e}  |  {rc.GetRate():>{wid}.{2}f} | {rc.mCorrBkgRate:>{wid}.{2}f}") for rc in self.mCoBare]

    return
  
  def buildOutArrays(self):

    self.gammaMap = {"csb": self.mCsBare, "cslp": self.mCsLP, "cob": self.mCoBare, "colp": self.mCoLP}
    return
  
  def findRunByHeight(self, height, key = "csb"):
    # [print(rc.mHeight) for rc in self.gammaMap.get(key)]
    found = next((rc for rc in self.gammaMap.get(key) if rc.mHeight == height), None)
    if found == None:
      print(f"ERROR: no mc object found in rc {key} with height {height}")
    return found
  
  
  # def plot_singleRows(self, fig: plt.figure, selection = "mCsBare",):
  #   Nrows = len(getattr(self, selection))

  #   for i, rc in enumerate(getattr(self, selection)):
  #     ax = fig.add_subplot.hist(Nrows, 1, i+1)
  #     ax.hist(x = rc.mBinX, bins = rc.mBin_edges, weights = rc.mBin_rawCounts, histtype = "step", label = "data " + str(hScale - rcBare_55cm.mHeight) + " cm")
  #   return
