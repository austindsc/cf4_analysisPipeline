#!/usr/local/bin/python3


"""
	this code directly plots SPE distributions from the data (good for checking SPE triggered data, with correct digitizer range)
"""


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
from ROOT import TMath, TFormula

print("Austin's makeshift waveform viewer")

# index of recorded waveforms in the event
iSipm = 0 # sipm data is first
iDis = 1 # discriminator data second

####### parsing #######
pars = argparse.ArgumentParser()
pars.add_argument("--first", type = int, required = False, default = 0, help = "what event to start in")
pars.add_argument("--g", type = int, required = False, default = 0, help = "what pulse group to start on")
pars.add_argument("--N", type = int, required = False, default = 0, help = "N events to plot together")
pars.add_argument("--live", type = int, required = False, default = False, help = "do viewing with user input")
pars.add_argument("--run", type = str, required = False, default = "20210726_0")
args=pars.parse_args()
live = bool(args.live)
run = args.run # LED data in CF4 gas?


# files and paths
raw_directory = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/rawdata"

######### BUILDING HISTS #########
qMax = 10
hq_range1 = TH1F("hq_range1", "hq_range1", 200, -1, qMax)
hq_smart1 = TH1F("hq_smart1", "hq_smart1", 200, -1, qMax)
hqc_smart1 = TH1F("hqc_smart1", "hqc_smart1", 200, -1, qMax)


# setup plotting
fig, axes = plt.subplots(1, 1, figsize = (8, 6))
ax0 = axes
if live == True:
	plt.ion()
	plt.show()

maxeV = 1 # how many event groups to look through

# loop through pulse groups (badly named event groupings, aka 'events')
for iev in range(args.g, maxeV):

	## load the data
	path = os.path.join(raw_directory, run, str(iev), "PMTtraces.bin")
	data = rb.ReadBlock(path)

	# pulse information
	traces = data.get('traces')
	vScale = data.get('v_scale')
	offset = data.get('v_offset')

	# loop through args.N pulses, or all if N = 0, (or less if not many events)
	for ind in range(args.first, args.first + (lambda: (len(traces), min(len(traces), args.N))[args.N > 0])()):
		
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

		if minTime < 300: # pulse is too early
			blines = getBaselineAvg(data = s1v, start = minTime+400, end = minTime+900)
			if live == True:
				print("using end as baseline")

		qRange, rEarly, rLate = minFixRangeIntegrate(data = s1v, bline = blines.get('v'), early = 150, late = 250)
		qSmart1, fEarly, fLate = toFlatRangeIntegrate(data  = s1v, bline = blines.get('v'), blineRMS = blines.get('e'), extra = 0)

		## try to do some fittinng ##

		## filling histograms
		hq_range1.Fill(qRange)
		hq_smart1.Fill(qSmart1)

		if live == True:
			#hold for user input
			print("\n event #", ind)
			print("minTime = ", minTime)
			print("blines = ", blines)
			print("qRange = ", qRange)
			print("qSmart = ", qSmart1)
			input()
			plt.pause(0.1)
			ax0.clear() # somehow reset but keep last N pulses

ax0.minorticks_on()
ax0.grid(True)
# ax0.set_ylim(-200, 200)

c1 = TCanvas( 'c1', 'c1', 600, 600)
hq_smart1.Draw("hist")

hq_range1.SetLineColor(3)
hq_range1.Draw("hist same")

c1.Modified()
c1.Update()

# plt.show()
ROOT.gApplication.Run()



