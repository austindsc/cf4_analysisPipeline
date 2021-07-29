#!/usr/local/bin/python3

import numpy as np
import os
import matplotlib.pyplot as plt
import argparse

import sys
sys.path.append('..')
from SBCcode.DataHandling import ReadBinary as rb

from ozTools import *

import ROOT
from ROOT import TFile
from ROOT import TCanvas, gApplication, gStyle, TLegend, TPaveText
from ROOT import TF1, TH1F, TGraph, TGraphErrors, TFitResult
from ROOT import TSpectrum
from ROOT import TMath, TFormula

# turn off all drawing stuff off 'live = 0'
# plot p0 


print("Austin's makeshift waveform viewer")

# index of recorded waveforms in the event
iSipm = 0 # sipm data is first
iDis = 1 # discriminator data second

pars = argparse.ArgumentParser()
pars.add_argument("--first", type = int, required = False, default = 0, help = "what event to start in")
pars.add_argument("--g", type = int, required = False, default = 0, help = "what pulse group to start on")
pars.add_argument("--N", type = int, required = False, default = 50, help = "N events to plot together")
pars.add_argument("--run", type = str, required = False, default = "20210726_0")
pars.add_argument("--live", type = int, required = False, default = False, help = "do viewing with user input")

args=pars.parse_args()
run = args.run # LED data in CF4 gas?
live = bool(args.live)

# files and paths
raw_directory = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/rawdata"

######### BUILDING HISTS #########
qMax = 10
hq_fAll = TH1F("hq_fAll", "hq_fAll", 200, -1, qMax)
hq_fChi2 = TH1F("hq_fChi2", "hq_fChi2", 200, -1, qMax)

hq_range1 = TH1F("hq_range1", "hq_range1", 200, -1, qMax)
hq_smart1 = TH1F("hq_smart1", "hq_smart1", 200, -1, qMax)

# setup plotting
fig, axes = plt.subplots(1, 1, figsize = (8, 6))
ax0 = axes

fitOpt = "WQN0"
specOpt = "goff"
if live == True:
	c1 = TCanvas('c1', 'c1', 600, 600)
	plt.ion()
	fitOpt = "WN0"
	specOpt = ""
	# plt.show()

maxeV = 1 # how many event groups to look through

# don't loop through pulse groups (badly named event groupings, aka 'events')
## load the data
path = os.path.join(raw_directory, run, str(args.g), "PMTtraces.bin")
data = rb.ReadBlock(path)

# pulse information
traces = data.get('traces')
vScale = data.get('v_scale')
offset = data.get('v_offset')

# loop through 1000 pulses (or less if not many events)
for ind in range(args.first, args.first + min(len(traces), args.N)):
	
	if ind%1000 == 0: # printing helper
		print("event", ind)

	# recorded sipm waveform
	s1 = traces[ind, iSipm][0:-100] # shape
	s1v = s1 * vScale[ind, iSipm] +  offset[ind, iSipm] # shape x scale + offset
	d1 = traces[ind, iDis] # discriminator

	# # trigger from discriminator
	# trigger = traces['traces'][ind,1] * traces['v_scale'][ind,1] + traces['v_offset'][ind,1]

	gap = 1 # pointless to multiply by dt, which is just 1ns
	# gap = data.get('dt')[ind, iSipm] * 10**9
	time = np.array(range(len(s1))) * gap #+ data.get("t0")[ind, iSipm]

	# ax0.plot(time, s1v, alpha = 0.4)
	# ax0.plot(time, d1)
	# plt.plot(sipm[:-100])
	# plt.plot(trigger[:-100]*0.1)	

	######## do processing ########
	minTime = np.argmin(s1v)
	blines = getBaselineAvg(data = s1v, start = 0, end = 200)

	spread = max(s1v) - min(s1v)
	print("spread = ", spread)
	# if spread < 0.007:
	# 	continue

	if minTime < 300: # pulse is too early
		blines = getBaselineAvg(data = s1v, start = minTime+400, end = minTime+900)
		if live == True:
			print("using end as baseline")

	qRange, rEarly, rLate = minFixRangeIntegrate(data = s1v, bline = blines.get('v'), early = 150, late = 250)
	qSmart1, fEarly, fLate = toFlatRangeIntegrate(data  = s1v, bline = blines.get('v'), blineRMS = blines.get('e'), extra = 0)

	# building TH1 (and flipping)
	hwf = TH1F("hwf", "hwf", len(s1v), time[0], time[-1])
	for i in range(len(s1v)):
		hwf.SetBinContent(i+1, -s1v[i])

	if live == True:
		# root drawing first so spectrum shows up
		hwf.Draw("hist")

	# baseline and peak finding
	spec = TSpectrum(1)
	ROOT.gVerboseLevel = ROOT.kError
	ROOT.gErrorIgnoreLevel = ROOT.kFatal
	spec.Search(hwf, 10., specOpt + "noMarkov", 0.5)

	###### fitting functions ######
	flat1 = TF1("flat1", "pol0", 0, 250)
	hwf.Fit(flat1, fitOpt, "goff")

	# 0 - height
	# 1 - rise time
	# 2 - fall time (Tl)
	# 3 - shift
	# 4 - baseline
	pulse1 = TF1("pulse1", "[0]/[2] * (TMath::Exp(-(x-[3])/([2] + [1])) - TMath::Exp(-(x - [3])/[1])) + [4]", spec.GetPositionX()[0] - 300, spec.GetPositionX()[0] + 800)
	pulse1.SetParameters(spec.GetPositionY()[0]*200, 70, 80, minTime, flat1.GetParameter(0))
	pulse1.SetParLimits(0, 0.5, 100)
	pulse1.SetParLimits(1, 10, 100)
	pulse1.SetParLimits(2, 40, 250)
	pulse1.SetParLimits(3, minTime - 150, minTime + 100)
	# pulse1.SetParLimits(4, -blines.get('v'), -blines.get('v'))
	pulse1.SetParLimits(4, flat1.GetParameter(0), flat1.GetParameter(0))

	hwf.Fit(pulse1, fitOpt, "goff", spec.GetPositionX()[0] - 100, spec.GetPositionX()[0] + 600)

	qFit = pulse1.GetParameter(0)
	chi2Fit = pulse1.GetChisquare()

	hq_fAll.Fill(qFit)
	if chi2Fit <= 2.0:
		hq_fChi2.Fill(qFit)

	## filling simple histograms
	hq_range1.Fill(qRange)
	hq_smart1.Fill(qSmart1)

	if live == True:
		# root drawing
		# hwf.Draw("hist")
		pulse1.Draw("same")
		flat1.Draw("hist same")

		c1.Modified()
		c1.Update()

		#hold for user input
		print("\n event #", ind)
		print("minTime = ", minTime)
		print("peaks found = ", spec.GetPositionX()[0])
		print("blines = ", blines)
		print("qRange = ", qRange)
		print("qSmart = ", qSmart1)
		print("Qfit = ", qFit, "    ch2 = ", chi2Fit, "   NDF = ", pulse1.GetNDF())
		input()
		# plt.show()
		# ax0.clear() # somehow reset but keep last N pulses

	del hwf

if live == False:
	c1 = TCanvas('c1', 'c1', 600, 600)

	hq_smart1.SetLineColor(1) # black
	hq_smart1.Draw("hist")

	hq_range1.SetLineColor(3) # green
	hq_range1.Draw("hist same")

	hq_fAll.SetLineColor(4) # blue
	hq_fAll.Draw("hist same")

	hq_fChi2.SetLineColor(6) # pink
	hq_fChi2.Draw("hist same")

	c1.Update()
	ROOT.gApplication.Run()



