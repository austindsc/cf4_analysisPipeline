#!/usr/local/bin/python3

#stock python stuff
import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
import time as chrono

## this is pyroot
from ROOT import ROOT
from ROOT import TVectorD
from ROOT import TCanvas, gApplication, TGraph, TPad, TRatioPlot, TLine, TText, TLatex, gStyle, TLegend
from ROOT import TGraphErrors, TF1, TFile, TH2F
from ROOT import gROOT

import sys
# sys.path.append('..')
from SBCcode.DataHandling import ReadBinary as rb

## my stuff ##
from ozTools import *

#########################
print("yoooo doing some SPE analysis")

#########################
# improvements from previous analysis code

# increase histogram range to 400
# store start and end time (seconds)


######### PARSING #########
pars = argparse.ArgumentParser()
pars.add_argument("--view", type = int, required = False, default = False)
pars.add_argument("--skip", type = int, required = False, default = 1)

pars.add_argument("--max", type = int, required = False, default = 1000000)
pars.add_argument("--outName", type = str, required = False)
pars.add_argument("--date", type = str, required = False)
pars.add_argument("--V", type = int, required = False)
pars.add_argument("--runFile", type = str, required = False)
pars.add_argument("--source", type = str, required = False)
args = pars.parse_args()

maxData = args.max
globalCount = 0 #  how many total events processed
skip = args.skip # how many  to skip through for pulse viewing

view = bool(args.view) # true (1) for viewing pulses

######### RUN LISTS, DATA INDICES, DIRECTORIES #########
path = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/"

outPath = path + "proc2/"
outName = path + "testSane2.root"
voltage = 53
ddate = "20201123"
runFile = "gasBig.txt"
source = 'bkg'

if args.outName != None:
	outName = str(outPath + args.outName + ".root")

if args.V != None:
	voltage = args.V

if args.date != None:
	ddate = str("2020" + args.date)

if args.runFile != None:
	runFile = args.runFile

if args.source != None:
	source = args.source

# items not in args parser
ampGain = 20
notes = "notes"

runList = np.loadtxt(path + "runLists/" + runFile, delimiter = ', ', dtype = str, unpack = False, skiprows = 1)

# j is the run txt file column
jRun = 0
jEv = 1
jSrc = 2
jV = 3
jP = 4
jT = 5
jPb = 6
jGain = 7
jOther = 8

# print(runList[5,0]) #(row, column), so runList[0, jT] gives temp of first run

# files and paths
raw_directory = path + "SBC-20-data/"

iSipm = 0 # sipm data is first
iDis = 1 # discriminator data second

############ setup plotting ############
if view == True:
	fig, axes = plt.subplots(1, 1, figsize = (8, 5))
	ax0 = axes
	plt.ion()
	plt.show()

# # root files
hFile = gROOT.FindObject(outName)
if hFile:
	hFile.Close()
hFile = TFile(outName, 'RECREATE', 'bigSane output file')


######### BUILDING HISTS #########
qFMax = 800
h2_hw = TH2F("h2_hw", "h2_hw", 2000, -0.01, 0.9, 2000, 0, 2000)
h2_hw_cuts = TH2F("h2_hw_cuts", "h2_hw_cuts", 2000, -0.01, 0.9, 2000, 0, 2000)

h2_qw = TH2F("h2_qw", "h2_qw", 8000, -10, qFMax, 2000, 0, 2000)
h2_qw_cuts = TH2F("h2_qw_cuts", "h2_qw_cuts", 8000, -10, qFMax, 2000, 0, 2000)

h2_hq = TH2F("h2_hq", "h2_hq", 2000, -0.01, 0.9, 8000, -10, qFMax)
h2_hq_cuts = TH2F("h2_hq_cuts", "h2_hq_cuts", 2000, -0.01, 0.9, 8000, -10, qFMax)

hF_all = TH1F('hF_all', 'hF_all', 8000, -10, qFMax)
hF_cuts = TH1F('hF_cuts', 'hF_cuts', 8000, -10, qFMax) # fixed range integral

hH_all = TH1F('hH_all', 'hH_all', 2000, -0.01, 0.9)
hH_cuts = TH1F('hH_cuts', 'hH_cuts', 2000, -0.01, 0.9)

h2_gq_raw = TH2F("h2_gq", "charge and time between all events", 10000, 0., 1., 8000, -10, qFMax)
h2_gq_slant = TH2F("h2_gq_slant", "charge and time between events, no slants", 10000, 0., 1., 8000, -10, qFMax)
h2_gq_cuts = TH2F("h2_gq_cuts", "charge and time between events, all cuts", 10000, 0., 1., 8000, -10, qFMax)

liveTime = 0.

######### looping on data #########
for line in runList: # getting all data in runlist

	if globalCount >= maxData:
		break

	### CUTS: env/daq level, voltage, source, etc
	if int(line[jV]) != voltage:
		continue
	if str(line[jSrc]) != source:
		continue
	if str(line[jRun][:8]) != ddate:
		continue
	if int(line[jGain]) != ampGain:
		continue
	if str(line[jOther]) != notes:
		continue

	# getting the path
	path = os.path.join(raw_directory, line[jRun], line[jEv], "PMTtraces.bin")
	data = rb.ReadBlock(path)
	print(path)

	# pulse information
	traces = data.get('traces')
	vScale = data.get('v_scale')
	offset = data.get('v_offset')

	# timing information
	tStart = data.get("t0_sec")[0, iSipm] + data.get("t0_frac")[0, iSipm]
	tEnd = data.get("t0_sec")[-1, iSipm] + data.get("t0_frac")[-1, iSipm]
	
	liveTime = liveTime + (tEnd - tStart)
	times = data.get("t0_sec")[:, iSipm] + data.get("t0_frac")[:, iSipm]

	# loop through 100 pulses of one trace object
	# ind loops through pulses
	print("N events =", len(traces), "\n")
	for ind in range(0, len(traces), skip):

		# time dif between pulses
		tGap = times[ind] - times[ind -1]

		globalCount += 1
		if globalCount >= maxData:
			break

		if ind%1000 == 0 and ind > 0:
			print("ev#", ind)

		if view == True:
			print("\n event #", ind)

		# if view == True:
		# 	print("t0_sec =", data.get("t0_sec")[ind, iSipm])
		# 	print("t0_frac =", data.get("t0_frac")[ind, iSipm])
		# 	print("t0 =", data.get("t0")[ind, iSipm])

		######### Getting Data #########
		# recorded sipm waveform
		s1 = traces[ind, iSipm][0:-50] # shape, last bits are wrong sometimes so skip last 30
		s1v = s1 * vScale[ind, iSipm] +  offset[ind, iSipm] # shape x scale + offset (actual pulse data)
		sMin = min(s1v)
		sMax = max(s1v)
		d1 = traces[ind, iDis][0:-50]*vScale[ind, iDis] + offset[ind, iDis] # discriminator
		gap = 1 # pointless to multiply by dt, which is just 1ns
		# gap = data.get('dt')[ind, iSipm] * 10**9
		time = np.array(range(len(s1v))) * gap # + data.get("t0")[ind, iSipm] # ignore the small shift
		minTime = np.argmin(s1v)

		######### processing #########

		# simple blines calc (FIRST 100 pts)
		# blines = getBaselineAvg(data = s1v, start = 0, end = 500)
		blines = getBaselineAvg(data = s1v, start = 0, end = 200)
		
		# recalculate baseline is pulse is early
		if minTime < 300:
			blines = getBaselineAvg(data = s1v, start = minTime+400, end = minTime+900)
			if view == True:
				print("using end as baseline")

		qRange, rEarly, rLate = minFixRangeIntegrate(data = s1v, bline = blines.get('v'), early = 250, late = 450)
		qFlat, fEarly, fLate = toFlatRangeIntegrate(data  = s1v, bline = blines.get('v'), blineRMS = blines.get('e'), extra = 0)
		qFlat2, fEarly2, fLate2 = toFlatRangeIntegrate(data  = s1v, bline = blines.get('v'), blineRMS = blines.get('e'), extra = 50) # good one
		#less than 1% error on large pulses from extra 50 ns of overshoot, but behaves much better for small pulses
		width = np.abs(fLate2 - fEarly2) - 100

		height = minAverage(data = s1v, bline = blines.get('v'), spread = 8)

		######### CUTTING and FILLING #########
		h2_gq_raw.Fill(tGap, qFlat2)

		# these pulses are just a flat line starting above 0 with a small area
		if sMax > 0 and sMin > -0.01 or (blines.get('v') - sMin) < 0.006: #slant cut (recovery pulses, have something above zero with no dip)
			if view == False:
				continue
			else:
				print("SLANT")

		h2_gq_slant.Fill(tGap, qFlat2)

		hF_all.Fill(qFlat2)
		hH_all.Fill(height)

		h2_hq.Fill(height, qFlat2)
		h2_qw.Fill(qFlat2, width)
		h2_hw.Fill(height, width)

		# 50 ns ahead of pulse return to baseline is at zero (start of pulse is cut off)
		if fEarly2 == 0:
			if view ==  True:
				print("actually too early")
			else:
				continue

		# bad small pulse (compare toFlat and fixed integral)
		if np.abs(qFlat2-qRange)/qFlat2 > 0.2 and (qFlat2 < 4 and qFlat2 > 0.5):
			if view == True:
				print("bad integral")
			else:
				continue

		h2_gq_cuts.Fill(tGap, qFlat2)

		hF_cuts.Fill(qFlat2)
		hH_cuts.Fill(height)

		h2_hq_cuts.Fill(height, qFlat2)
		h2_qw_cuts.Fill(qFlat2, width)
		h2_hw_cuts.Fill(height, width)

		######### Viewing/Troubleshooting #########
		if view == True:

			#  only view certain size waveforms
			# if qFlat2 < 20 or qFlat2 > 60:
			# 	continue

			# for looking at pulses
			ax0.plot(time, s1v, '.', alpha = 0.7)
			ax0.plot(time, d1, '-') # discriminator pulse
			ax0.axvline(rEarly, color = 'k', label = "fixed Range")
			ax0.axvline(rLate, color = 'k')

			ax0.axvline(fEarly, color = 'r', label = "variable Range")
			ax0.axvline(fLate, color = 'r')

			ax0.axvline(fEarly2, color = 'g', label = "variable2 Range")
			ax0.axvline(fLate2, color = 'g')

			ax0.plot(time, blines.get('v')*np.ones(len(time)), label = "early Baseline")
			# ax0.plot(time, blinesBack.get('v')*np.ones(len(time)))
			ax0.legend()

			print("ev#", ind)
			# print(blines.get('v'), blines.get('v') - sMin)
			# print( '%0.1E' % blines['v'], '%0.1E' % blines['e'])
			# print('%0.1E' % blinesH[0], '%0.1E' % blinesH[1])
			print("%i width" % width)
			print("%0.2f fixed q" % qRange)
			print("%0.2f var q" % qFlat)
			print("%0.2f var2 q" %qFlat2, "\n")

			input()
			# plt.pause(0.4)
			ax0.clear()

hFile.cd()

rootLiveTime = TVectorD(1)
rootLiveTime[0] = liveTime
rootLiveTime.Write("liveTime")

h2_gq_raw.Write()
h2_gq_slant.Write()
h2_gq_cuts.Write()

hF_all.Write()
hH_all.Write()

h2_hq.Write()
h2_qw.Write()
h2_hw.Write()

hF_cuts.Write()
hH_cuts.Write()

h2_hq_cuts.Write()
h2_qw_cuts.Write()
h2_hw_cuts.Write()

hFile.Close()