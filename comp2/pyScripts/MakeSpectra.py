"""
MakeSpectra - simple bookkeeping for CF4 spectral analysis.

Usage: python MakeSpectra.py fileConfig.json anaConfig.json

you can pass a single runID or a .txt file with list of runIDs to fileConfig["runList"]

"""

"""
Helpers copied over from analysis1
"""

import os
from re import VERBOSE
import sys
import numpy as np

try:
  import json5 as json
except ImportError:
  import json
  
import ROOT
# ROOT.gROOT.SetBatch(True) # to avoid annoying NSA app warning
ROOT.gROOT.SetBatch(False)
# ROOT.enableJSVis()
ROOT.gErrorIgnoreLevel = ROOT.kError # suppress warnings for peak buffer

## Atown's stuff
sys.path.append("../") # to get globalsC2
from globalsC2 import *

import helpers.anaStructs as ANAS # includes SORTING class/namespace
import helpers.getSPE as GSPE
import helpers.graphBuilders as GB

class AnalysisConfig:
  """Class to hold analysis configuration parameters."""
  def __init__(self, config_dict):
    strAna = "analysis"
    
    self.inputDict = config_dict.get(strAna)
    # Data analysis parameters
    self.nChan = config_dict.get(strAna).get("nchan")
    self.verbose = config_dict.get(strAna).get("verbose")
    self.qvar = config_dict.get(strAna).get("qvar")
    self.peCut = config_dict.get(strAna).get("peCut")
    self.speWidth = config_dict.get(strAna).get("speWidth")
    
    # Global analysis parameters
    self.binLow = config_dict.get(strAna).get("binLow")
    self.binHigh = config_dict.get(strAna).get("binHigh")
    self.nBins = config_dict.get(strAna).get("nBins")
    self.binWidth = (self.binHigh - self.binLow) / self.nBins if self.nBins else None
    
    self.useRate = config_dict.get(strAna).get("useRate", True)
    self.doAverage = config_dict.get(strAna).get("sumRuns", False)
    
    # for other higher level analysis
    # strGlobal = "global_analysis"
    pass

class FileConfig:
  """Class to hold file configuration parameters."""
  def __init__(self, config_dict):
    strFiles = "files"
    self.inputDict = config_dict.get(strFiles)
    self.dataRoot = config_dict.get(strFiles).get("dataRoot")
    self.dataRoot = self.dataRoot.lstrip('/')
    
    self.runList = config_dict.get(strFiles).get("runList")
    self.filePattern = config_dict.get(strFiles).get("filePattern", "npv1_{}.root")
    self.treeName = config_dict.get(strFiles).get("treeName", "tt")
    pass


def th1_to_numpy(h):
  """Safe conversion of a ROOT TH1 to numpy arrays.
  Returns (edges, counts, errors). Avoids ctypes buffer views to prevent segfaults.
  Uses variable bin edges if present, else per-bin low edges.
  """
  nb = h.GetNbinsX()
  # Bin edges: prefer explicit bin low edges to avoid axis bins quirks
  edges = np.array([h.GetBinLowEdge(i + 1) for i in range(nb + 1)], dtype=float)
  # Counts and errors using vectorized fromiter (safe with PyROOT)
  counts = np.fromiter((h.GetBinContent(i + 1) for i in range(nb)), dtype=float, count=nb)
  errors = np.fromiter((h.GetBinError(i + 1) for i in range(nb)), dtype=float, count=nb)
  return edges, counts, errors

## 
def ExtractSPE(fileConfig, anaConfig='analysis.json', argRunList = None):
  """Extract SPE values for all runs without building full histograms.
  
  Returns:
    runIDs: list of run IDs
    speValues: list of [spe_ch0, spe_ch1, spe_ch2] for each run
  """
  # Load configs
  with open(fileConfig) as f:
    fileCfgDict = json.load(f)
  with open(anaConfig) as f:
    anaDict = json.load(f)

  anaCfg = AnalysisConfig(anaDict)
  fileCfg = FileConfig(fileCfgDict)
  print(fileCfg.inputDict)

  # Resolve data path and run list
  dataPath = fileCfg.dataRoot
  file_pattern = fileCfg.filePattern
  runs = fileCfg.inputDict.get("runList", {})
  if argRunList is not None:
    runs = argRunList
  elif type(runs) is str:
    runs = np.loadtxt(os.path.join(COMP2PATH, "runLists", fileCfg.runList), dtype=str).tolist()
  elif type(runs) is list:
    runs = runs
  elif type(runs) is int:
    runs = [runs]
  print("getting runs:", runs) if VERBOSE else None
  
  # File checks
  existing = []
  MakePath = lambda run: os.path.join(CF4HOME, dataPath, file_pattern.format(run))
  for run in runs:
    if os.path.isfile(MakePath(run)):
      existing.append(run)

  # Extract SPE for each run
  runIDs = []
  speValues = []
  tFiles = []
  rcObjList = []

  for run in existing:
    fpath = MakePath(run)
    f = ROOT.TFile.Open(fpath)
    if not f or f.IsZombie():
      continue
    tFiles.append(f)
    tt = f.Get(fileCfg.treeName)
    if tt is None:
      continue

    rc = ANAS.RunClass(run, tt, fname=fpath)
    rc.computeSPE(width=anaCfg.speWidth, verbose=anaCfg.verbose, varname=anaCfg.qvar, channels=[i for i in range(anaCfg.nChan)], nBins = anaCfg.nBins)
    
    rcObjList.append(rc)
    runIDs.append(run)
    speValues.append([rc.GetChan(i).mSpe for i in range(anaCfg.nChan)])

  return runIDs, speValues, rcObjList, tFiles


# def CallUprootSetup(rcList):
#   """Call UprootConstructor for each RunClass in rcList."""
#   for rc in rcList:
#     rc.UprootConstructor()

def Build_histogram(fileConfig, anaConfig='analysis.json'):
  """Light refactor: load configs, check files, construct a ROOT TH1 and a numpy container.

  Returns (hReference, hNp) where:
    - hReference: ROOT.TH1F configured using global analysis bins
    - hNp: dict with keys {bin_edges, counts, errors}
  """
  # Load configs
  with open(fileConfig) as f:
    fileCfgDict = json.load(f)
  with open(anaConfig) as f:
    anaDict = json.load(f)

  # Analysis parameters (use values directly from AnalysisConfig)
  anaCfg = AnalysisConfig(anaDict)
  fileCfg = FileConfig(fileCfgDict)

  # Resolve data path and run list using FileConfig class
  dataPath = fileCfg.dataRoot
  file_pattern = fileCfg.filePattern
  runs = fileCfg.inputDict.get("runList")
  # check if runs isnt a list
  if not isinstance(runs, list):
    runs = [runs]
  if isinstance(runs, str):
    runs = np.loadtxt(os.path.join(COMP2PATH, "runLists", fileCfg.runList), dtype=str).tolist()

  # File checks (optional bookkeeping)
  existing, missing = [], []
  MakePath = lambda run: os.path.join(CF4HOME, dataPath, file_pattern.format(run))
  for run in runs:
    path = MakePath(run)
    relPath = os.path.join(dataPath, file_pattern.format(run))
    if os.path.isfile(path):
      existing.append(run)
    else:
      # Keep lightweight—report missing but continue
      print("Missing file:", relPath)
      missing.append(run)

  # Construct ROOT histogram (empty here; downstream code may fill)
  hReference = ROOT.TH1F("hReference", "hReference", anaCfg.nBins, anaCfg.binLow, anaCfg.binHigh)
  hReference.Sumw2(True)

  # === Core loop: open files, build RunClass, compute SPE, and fill spectrum ===
  tFiles = []  # keep ROOT files alive
  rcList = []
  speList = []
  
  for run in existing:
    fpath = MakePath(run)
    f = ROOT.TFile.Open(fpath)
    if not f or f.IsZombie():
      print(f"Warning: failed to open {fpath}")
      continue
    tFiles.append(f)
    tt = f.Get(fileCfg.treeName)
    if tt is None:
      print(f"Warning: missing tree '{fileCfg.treeName}' in {fpath}")
      continue
    ## actual physics
    rc = ANAS.RunClass(run, tt)
    rc.computeSPE(width=anaCfg.speWidth, verbose=anaCfg.verbose, varname=anaCfg.qvar, nBins = 200, channels=[i for i in range(anaCfg.nChan)])
    rc.UprootConstructor()
    rcList.append(rc)
    speList.append([[chan.mSpe.spe, chan.mSpe.chi2] for chan in rc.mChannels]) # for possible bookkeeping
    
    rc.BuildSpectrum(histClone=hReference, varname=anaCfg.qvar, peCut=anaCfg.peCut, rate=anaCfg.useRate)

  if anaCfg.verbose:
    input("say somthing to move past bkg SPE plots")

  ## do averaging if needed
  hAvg = None
  if anaCfg.doAverage and len(rcList) > 1:
    print("doing averaging over", len(rcList), "runs")
    hAvg = hReference.Clone("hAvg")
    GB.build_bkgHistogramFast(rcList, anaCfg.qvar, hAvg, peCut=anaCfg.peCut, verbose=anaCfg.verbose)
  
  ## select the correct return histogram (robust)
  if len(rcList) == 1:
    hReturn = rcList[0].mhSpectrum.Clone("hOut")
  else:
    hReturn = hAvg.Clone("hOut")
  hReturn.SetDirectory(0) # dangling objects are dangerous in root haha
  
  print("Done building histogram:", hReturn.GetName())
  
  # making numpy histogram from ROOT one (fast, vectorized)
  binEdges, counts, errors = th1_to_numpy(hReturn)
  hNp = {'bin_edges': binEdges, 'counts': counts, 'errors': errors}
  return hReturn, hNp


def main(configFile, anaConfig='analysis.json'):
  """Main function to build histogram from config files."""
  # Call the refactored function so script remains runnable standalone
  hOut, hNp = Build_histogram(configFile, anaConfig)
  # Minimal confirmation output

  ## np printout
  print("numpy histogram output!")
  print("TH1F type:", type(hOut))
  
if __name__ == '__main__':
  if len(sys.argv) < 2:
    print("Usage: python MakeSpectra.py fileConfig.json anaConfig.json")
    sys.exit(1)
  main(sys.argv[1], sys.argv[2])