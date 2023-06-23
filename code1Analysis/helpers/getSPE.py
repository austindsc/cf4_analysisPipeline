
import ROOT
import numpy as np
from array import array
# from helpers.anaStructs import SpeStuff

## ---------------------------------------------------------------- ##
class SpeStuff:

  def __init__(self, spe = -1, speE = -1, chi2 = 0, rawHist = None):
    self.spe = spe
    self.espe = speE
    self.chi2 = chi2 # maximum value
    self.mRawHist = rawHist
    return

  def print(self):
    print(f"spe = {self.spe:0.2f} +/- {self.espe:0.2f}     maxChi2 = {self.chi2:0.2f}")
    return

def getSPEVval2(tt, width = 10., verbose = False, varName = "qFixed", channel = 0, runName = "", nBins = 100, specialCase = 0):

  histName = "qvar-hist_"+runName+"_chan"+str(int(channel))
  drawCommand = str(varName+f"[{channel:.0f}]" + ">>"+histName)
  cutCommand = str(varName+f"[{channel:.0f}]"+"> 30")
  dopt = "goff"
  fopt = "QRN0ML"
  if verbose == True:
    dopt = ""
    fopt = "+RML same"
    c = ROOT.TCanvas(drawCommand, drawCommand)
    c.cd()
  
  findPeakFactor = 1.4

  hist = ROOT.TH1F(histName, histName, nBins, 0, 1000)
  hist.GetXaxis().SetTitle("Charge, adc x ns")
  hist.GetYaxis().SetTitle("Counts")
  # hist.GetXaxis().SetRange(int(100/5), 200)

  tt.Draw(drawCommand, cutCommand, dopt)

  ## First peak as maximum
  xMaxPeak = hist.GetXaxis().GetBinCenter(hist.GetMaximumBin())
  if verbose:
    print("max location = ", xMaxPeak)
  fg1 = ROOT.TF1("fg1", "gaus", xMaxPeak - width, xMaxPeak + width)
  fg1.SetParameters(hist.GetMaximum(), xMaxPeak, width/2)
  hist.Fit("fg1", fopt)
  sig1pe = fg1.GetParameter(2)
  q1pe = fg1.GetParameter(1)

  ## second peak
  xa = hist.GetXaxis()
  xa.SetRange(xa.FindBin(2*xMaxPeak - 4*width), xa.FindBin(2*xMaxPeak + 4*width))
  npeaks = 2
  sigma = 2.
  s = ROOT.TSpectrum(npeaks)
  s.Search(hist, sigma, "nobackground", 0.001)

  # s.Print()
  xp = [s.GetPositionX()[i] for i in range(npeaks)]
  if verbose:
    print("second peak = ", xp[0])
  xa.SetRange() # reset range
  peak2 = xp[0]
  if peak2 < 1.5*q1pe:
    peak2 = xp[1]

  if specialCase == 202204141154: # the 52V run
    xp = sorted(xp)
    peak2 = xp[0] # use second peak (first along x)

  ## fit second gaussian
  fg2 = ROOT.TF1("fg2", "gaus", peak2 - width, peak2 + width)
  heightGuess = hist.GetBinContent(hist.GetXaxis().FindBin(peak2))
  fg2.SetParameters(heightGuess, peak2, sig1pe)
  hist.Fit("fg2", fopt)
  q2pe = fg2.GetParameter(1)

  spe = q2pe - q1pe
  eSPE = np.sqrt(fg2.GetParError(1)**2 + fg1.GetParError(1)**2)
  speObject = SpeStuff(spe, eSPE, chi2 = max(fg1.GetChisquare(), fg2.GetChisquare()), rawHist = hist)

  if verbose == True:
    print("q1pe = ", q1pe)
    print("q2pe = ", q2pe)
    print(f"ch{channel:.0f} SPE, {spe:0.2f}")

    # ROOT.gStyle.SetOptFit(1)
    c.SetLogy()
    c.Update()
  else:
    del fg1, fg2, hist

  return speObject

def getSPEVval(tt, width = 10., verbose = False, varName = "qFixed", channel = 0, runName = ""):

  drawCommand = str(varName+f"[{channel:.0f}]" + ">>hist")
  dopt = "goff"
  fopt = "QRN0ML"
  if verbose == True:
    dopt = ""
    fopt = "+RML same"
    c = ROOT.TCanvas(runName+"_"+drawCommand, runName+"_"+drawCommand)
    c.cd()
  
  findPeakFactor = 1.4

  hist = ROOT.TH1F("hist", "hist", 200, 0, 1000)
  hist.GetXaxis().SetTitle("Charge, adc x ns")
  hist.GetYaxis().SetTitle("Counts")
  hist.GetXaxis().SetRange(int(100/5), 200)

  tt.Draw(drawCommand, "", dopt)

  ## Find first peak as maximum ##
  npeaks = 2
  s = ROOT.TSpectrum(npeaks)
  s.Search(hist, 1., dopt, 0.01)

  # s.Print()
  # xpeaks = s.GetPositionX()
  # yp = s.GetPositionY()
  xp = [s.GetPositionX()[i] for i in range(npeaks)]
  # ypeaks = [i for s.GetPositionY()[i] in range(npeaks)]
  # xp.sort()

  # if there is a low PE pedestal
  print(xp[1] / xp[0])
  if xp[1] < xp[0]:
    xp[1] = 1.95*xp[0]

  # print(xp)
  # print(xpeaks)

  ## fit three gaussians
  fg1 = ROOT.TF1("fg1", "gaus", xp[0] - width, xp[0] + width)
  fg1.SetParameters(hist.GetMaximum(), xp[0], width/2)
  hist.Fit("fg1", fopt)
  q1pe = fg1.GetParameter(1)

  fg2 = ROOT.TF1("fg2", "gaus", xp[1] - width, xp[1] + width)
  fg2.SetParameters(hist.GetMaximum()/20, xp[1], width/2)
  hist.Fit("fg2", fopt)
  q2pe = fg2.GetParameter(1)

  spe = q2pe - q1pe
  eSPE = np.sqrt(fg2.GetParError(1)**2 + fg1.GetParError(1)**2)

  speObject = SpeStuff(spe, eSPE, chi2 = max(fg1.GetChisquare(), fg2.GetChisquare()))

  if verbose == True:
    print("q1pe = ", q1pe)
    print("q2pe = ", q2pe)
    print(f"ch{channel:.0f} SPE, {spe:0.2f}")

    ROOT.gStyle.SetOptFit(1)
    c.SetLogy()
    c.Update()
  else:
    del fg1, fg2, hist

  return speObject