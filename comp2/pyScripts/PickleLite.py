## this is used to open the pickled objects saved by DataAnalysis all spectra

# Helper: load sanitized snapshots and reconstruct lightweight spectra objects
import os
import pickle
from typing import Any, Dict

from helpers.anaStructs import RunClass
from globalsC2 import *

import numpy as np
import matplotlib.pyplot as plt

## for opening ##
def load_snap_list(path):
  with open(path, "rb") as f:
    data = pickle.load(f)
  return [LiteRun(snap) for snap in data]

# Example usage (uncomment and set COMP2PATH if not already in scope):
# save_dir = os.path.join(COMP2PATH, "DataAnalysis", "saved_spectra")
# co_bare = load_snap_list(os.path.join(save_dir, "rc_Co_bare_source.snap.pkl"))
# print(type(co_bare[0]))
# print(f"Loaded {len(co_bare)} Co bare LiteRun objects. First IDs:", [lr.id for lr in co_bare[:3]])


# Helper: convert RunClass -> pure-Python dict without any ROOT/CPyCppyy objects
def runclass_snapshot(rc: RunClass) -> Dict[str, Any]:
  snap: Dict[str, Any] = {}
  snap["id"] = rc.id
  snap["name"] = getattr(rc, "mName", None)
  snap["height"] = getattr(rc, "mHeight", None)
  snap["runtime"] = getattr(rc, "mRunTime", None)
  snap["mNEvents"] = getattr(rc, "mNEvents", None)
  
  # Binning
  snap["bin_edges"] = np.array(rc.mBin_edges) if hasattr(rc, "mBin_edges") else None
  snap["bin_centers"] = np.array(rc.mBinX) if hasattr(rc, "mBinX") else None
  
  ## spectra
  snap["counts_raw"] = np.array(rc.mBin_rawCounts) if hasattr(rc, "mBin_rawCounts") else None
  snap["rate_raw"] = np.array(rc.mBin_rate) if hasattr(rc, "mBin_rate") else None
  snap["erate_raw"] = np.array(rc.mBin_erate) if hasattr(rc, "mBin_erate") else None
  
  # Background-subtracted (if present)
  snap["BkgCorrected_rawCounts"] = np.array(rc.BkgCorrected_rawCounts) if hasattr(rc, "BkgCorrected_rawCounts") else None
  snap["BkgCorrected_errors"] = np.array(rc.BkgCorrected_errors) if hasattr(rc, "BkgCorrected_errors") else None
  
  # Channels SPE info
  try:
    chans = []
    anaCfg = MS.AnalysisConfig(os.path.join(COMP2PATH, "pyScripts/analysis.json"))
    for ich in range(anaCfg.nChan):
      ch = rc.GetChan(ich)
      chans.append({
        "index": ich,
        "spe": getattr(ch.mSpe, "spe", None),
        "espe": getattr(ch.mSpe, "espe", None),
        "gain": getattr(ch, "mGain", None),
      })
    snap["channels"] = chans
  except Exception:
    snap["channels"] = None
  return snap

##### for easy use and data manipulation in comparison with MC
class LiteRun:
  def __init__(self, snap: dict):
    self.id = snap.get("id")
    self.name = snap.get("name")
    # Match RunClass / RealDataSet naming (mHeight), keep .height for backwards compatibility
    self.mHeight = snap.get("height")
    self.height = self.mHeight
    self.channels = snap.get("channels")
    # RunClass uses GetRunTime(); we store raw value as runtime
    self.runtime = snap.get("runtime")
    self.mNEvents = snap.get("mNEvents")
    # Histogram binning: align names with RunClass (mBin_edges, mBinX)
    self.bin_edges = np.array(snap.get("bin_edges")) if snap.get("bin_edges") is not None else None
    self.bin_centers = np.array(snap.get("bin_centers")) if snap.get("bin_centers") is not None else None
    self.mBin_edges = self.bin_edges
    self.mBinX = self.bin_centers
    
    # Raw spectra: align with RunClass (mBin_rawCounts, mBin_rate, mBin_erate)
    self.counts_raw = np.array(snap.get("counts_raw")) if snap.get("counts_raw") is not None else None
    self.rate_raw = np.array(snap.get("rate_raw")) if snap.get("rate_raw") is not None else None
    self.erate_raw = np.array(snap.get("erate_raw")) if snap.get("erate_raw") is not None else None
    self.mBin_rawCounts = self.counts_raw
    self.mBin_rate = self.rate_raw
    self.mBin_erate = self.erate_raw
    
    self.BkgCorrected_rawCounts = np.array(snap.get("BkgCorrected_rawCounts")) if snap.get("BkgCorrected_rawCounts") is not None else None
    self.BkgCorrected_errors = np.array(snap.get("BkgCorrected_errors")) if snap.get("BkgCorrected_errors") is not None else None
    
    self.BkgCorrected_rate = self.BkgCorrected_rawCounts/self.runtime
    self.BkgCorrected_rateError = self.BkgCorrected_errors/self.runtime

  ## only background corrected plotting for now
  def mpl_errorbar(self, ax, useCounts=True, label=None, verbose = 0, **kwargs):
    if (useCounts and self.BkgCorrected_rawCounts is not None):
      y = self.BkgCorrected_rawCounts
      if verbose >= 1:
        print(f"LiteRun {self.id}: plotting background-corrected counts")
    else:
      y = self.BkgCorrected_rate
      if verbose >= 1:
        print(f"LiteRun {self.id}: plotting background-corrected rate")
    ey = self.BkgCorrected_errors if (useCounts and self.BkgCorrected_errors is not None) else self.BkgCorrected_rateError
    x = self.bin_centers
    if x is None or y is None:
      raise RuntimeError(f"LiteRun {self.id}: missing data to plot")
    a = ax.errorbar(x, y, yerr=ey, label=label or f"{self.height}mm", **kwargs)
    return a[0]  # return the Line2D object

  def mpl_hist(self, ax, useCounts=True, histtype='step', label=None, verbose=0, **kwargs):
    if (useCounts and self.BkgCorrected_rawCounts is not None):
      values = self.BkgCorrected_rawCounts
      if verbose >= 1:
        print(f"LiteRun {self.id}: plotting background-corrected counts")
    else:
      values = self.BkgCorrected_rate
      if verbose >= 1:
        print(f"LiteRun {self.id}: plotting background-corrected rate")
    
    if self.bin_edges is None or values is None:
      raise RuntimeError(f"LiteRun {self.id}: missing histogram data")
    n, bins, patches = ax.hist(self.bin_centers, bins=self.bin_edges, weights=values, histtype=histtype, label = label or f"{self.height}mm", **kwargs)
    return patches[0]
  
  def GetRate(self):
    return np.sum(self.rate_raw)
  
  def GetCorrectedRate(self):
    return np.sum(self.BkgCorrected_rate)
  
  def print(self):
    print(f"LiteRun id: {self.id}, name: {self.name}, height: {self.mHeight} cm, runtime: {self.runtime//60} m {int(self.runtime%60):d} s, NEvents: {self.mNEvents}")
    return
  
  def GetRunTime(self):
    return self.runtime

  def rebin(self, new_edges):
    """Rebin all histogram arrays to new bin edges.
    new_edges must span the same range or a subset of the original edges.
    Counts are summed; rates and errors are propagated appropriately.
    """
    old_edges = self.bin_edges
    if old_edges is None:
      raise RuntimeError(f"LiteRun {self.id}: no bin_edges to rebin from")

    new_centers = 0.5 * (new_edges[:-1] + new_edges[1:])
    nNew = len(new_centers)

    def _rebin_counts(old_counts):
      if old_counts is None:
        return None
      new_counts = np.zeros(nNew)
      for j in range(nNew):
        mask = (self.bin_centers >= new_edges[j]) & (self.bin_centers < new_edges[j+1])
        new_counts[j] = np.sum(old_counts[mask])
      return new_counts

    def _rebin_errors(old_errors):
      if old_errors is None:
        return None
      new_errors = np.zeros(nNew)
      for j in range(nNew):
        mask = (self.bin_centers >= new_edges[j]) & (self.bin_centers < new_edges[j+1])
        new_errors[j] = np.sqrt(np.sum(old_errors[mask]**2))
      return new_errors

    # Rebin all count-like arrays
    self.counts_raw = _rebin_counts(self.counts_raw)
    self.BkgCorrected_rawCounts = _rebin_counts(self.BkgCorrected_rawCounts)

    # Rebin errors in quadrature
    self.erate_raw = _rebin_errors(self.erate_raw)
    self.BkgCorrected_errors = _rebin_errors(self.BkgCorrected_errors)

    # Recompute rates from rebinned counts
    self.rate_raw = self.counts_raw / self.runtime if self.counts_raw is not None else None
    self.BkgCorrected_rate = self.BkgCorrected_rawCounts / self.runtime if self.BkgCorrected_rawCounts is not None else None
    self.BkgCorrected_rateError = self.BkgCorrected_errors / self.runtime if self.BkgCorrected_errors is not None else None

    # Update bin info
    self.bin_edges = new_edges
    self.bin_centers = new_centers
    self.mBin_edges = self.bin_edges
    self.mBinX = self.bin_centers
    self.mBin_rawCounts = self.counts_raw
    self.mBin_rate = self.rate_raw
    self.mBin_erate = self.erate_raw

  # def GetNP_integral(self, xLow, xHigh, option="rate"):
  #   """Get integral of numpy histogram between xLow and xHigh."""
  #   bin_mask = (self.bin_centers >= xLow) & (self.bin_centers <= xHigh)
  #   integral = np.sum(self.counts_raw[bin_mask]) if option == "counts" else np.sum(self.rate_raw[bin_mask])
  #   return integral