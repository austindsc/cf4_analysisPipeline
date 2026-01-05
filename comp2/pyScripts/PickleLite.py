## this is used to open the pickled objects saved by DataAnalysis all spectra

# Helper: load sanitized snapshots and reconstruct lightweight spectra objects
import os
import pickle
from typing import Any, Dict

from helpers.anaStructs import RunClass
from globalsC2 import *

import numpy as np
import matplotlib.pyplot as plt


# Helper: convert RunClass -> pure-Python dict without any ROOT/CPyCppyy objects
def runclass_snapshot(rc: RunClass) -> Dict[str, Any]:
  snap: Dict[str, Any] = {}
  snap["id"] = rc.id
  snap["name"] = getattr(rc, "mName", None)
  snap["height"] = getattr(rc, "mHeight", None)
  snap["runtime"] = getattr(rc, "mRunTime", None)
  snap["nEvents"] = getattr(rc, "mNEvents", None)
  # Binning + spectra
  snap["bin_edges"] = np.array(rc.mBin_edges) if hasattr(rc, "mBin_edges") else None
  snap["bin_centers"] = np.array(rc.mBinX) if hasattr(rc, "mBinX") else None
  snap["counts_raw"] = np.array(rc.mBin_rawCounts) if hasattr(rc, "mBin_rawCounts") else None
  snap["rate_raw"] = np.array(rc.mBin_rate) if hasattr(rc, "mBin_rate") else None
  snap["erate_raw"] = np.array(rc.mBin_erate) if hasattr(rc, "mBin_erate") else None
  # Background-subtracted (if present)
  snap["counts_bkgsub"] = np.array(rc.mBin_rawCounts_bkgCorr) if hasattr(rc, "mBin_rawCounts_bkgCorr") else None
  snap["rate_bkgsub"] = np.array(rc.mBin_rate_bkgCorr) if hasattr(rc, "mBin_rate_bkgCorr") else None
  snap["erate_bkgsub"] = np.array(rc.mBin_erate_bkgCorr) if hasattr(rc, "mBin_erate_bkgCorr") else None
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
    self.height = snap.get("height")
    self.runtime = snap.get("runtime")
    self.nEvents = snap.get("nEvents")
    self.bin_edges = np.array(snap.get("bin_edges")) if snap.get("bin_edges") is not None else None
    self.bin_centers = np.array(snap.get("bin_centers")) if snap.get("bin_centers") is not None else None
    self.counts_raw = np.array(snap.get("counts_raw")) if snap.get("counts_raw") is not None else None
    self.rate_raw = np.array(snap.get("rate_raw")) if snap.get("rate_raw") is not None else None
    self.erate_raw = np.array(snap.get("erate_raw")) if snap.get("erate_raw") is not None else None
    self.counts_bkgsub = np.array(snap.get("counts_bkgsub")) if snap.get("counts_bkgsub") is not None else None
    self.rate_bkgsub = np.array(snap.get("rate_bkgsub")) if snap.get("rate_bkgsub") is not None else None
    self.erate_bkgsub = np.array(snap.get("erate_bkgsub")) if snap.get("erate_bkgsub") is not None else None
    self.channels = snap.get("channels")

  def mpl_errorbar(self, ax, use_bkgsub=True, label=None, marker="o", linestyle="None"):
    y = self.rate_bkgsub if (use_bkgsub and self.rate_bkgsub is not None) else self.rate_raw
    ey = self.erate_bkgsub if (use_bkgsub and self.erate_bkgsub is not None) else self.erate_raw
    x = self.bin_centers
    if x is None or y is None:
      raise RuntimeError(f"LiteRun {self.id}: missing data to plot")
    ax.errorbar(x, y, yerr=ey, fmt=marker, linestyle=linestyle, label=label or f"{self.height}mm")

  def mpl_hist(self, ax, use_bkgsub=False, histtype='step', color=None, label=None):
    values = self.rate_bkgsub if (use_bkgsub and self.rate_bkgsub is not None) else self.rate_raw
    if self.bin_edges is None or values is None:
      raise RuntimeError(f"LiteRun {self.id}: missing histogram data")
    ax.hist(self.bin_centers, bins=self.bin_edges, weights=values, histtype=histtype, color=color, label=label)

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