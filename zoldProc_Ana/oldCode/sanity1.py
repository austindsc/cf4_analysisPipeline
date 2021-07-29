#!/usr/local/bin/python3

#stock python stuff
import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
import time as chrono

## this is pyroot
from ROOT import ROOT
from ROOT import TCanvas, TGraph, TPad, TRatioPlot, TLine, TText, TLatex, gStyle, TLegend
from ROOT import TGraphErrors, TF1, TFile
from ROOT import gROOT

import sys
# sys.path.append('..')
from SBCcode.DataHandling import ReadBinary as rb

## my stuff ##
from ozTools import *

#########################
print("yoooo doing some SPE analysis")
# 1. load in all runs (date and number)
# 2. get N events? or just try until read fails
# 3. search 

# can I remove noise from pulses by looking at noise in the disc. signal??

#### plans? ####
# okay lets look at SPE distributions for small runs and big runs

#compare low PE spectra to get better QPE values

# goal: find what 1PE is
# maybe we need to apply a moving average/filter
#########################

pars = argparse.ArgumentParser()
pars.add_argument("--view", type = int, required = False, default = False)
pars.add_argument("--max", type = int, required = False, default = 100000)
args = pars.parse_args()

maxData = args.max
view = bool(args.view) # true (1) for viewing pulses

## Getting run lists ##
runPath = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/atHome/runLists/"
runFile = "vacBig.txt"

voltage = 55
source = 'bkg'

runList = np.loadtxt(runPath + runFile, delimiter = ', ', dtype = str, unpack = False, skiprows = 1)

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
raw_directory = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/atHome/SBC-20-data/"

iSipm = 0 # sipm data is first
iDis = 1 # discriminator data second

############ setup plotting ############
if view == True:
	fig, axes = plt.subplots(1, 1, figsize = (8, 5))
	ax0 = axes
	plt.ion()
	plt.show()

# root files
hFile = gROOT.FindObject('SPE1.root')
if hFile:
	hFile.Close()
hFile = TFile('SPE1.root', 'RECREATE', 'speAna1 test file')

# Create a new canvas, and customize it.
# hq_noSlant = TH1F( 'hq_noSlant', 'sipm Q, bline Slant ', 400, -10, 20) # full integral
# hq_all = TH1F('hq_all', 'hq_all', 400, -10, 20)
c1 = TCanvas("c1", "c1", 600, 600)
c1.Divide(1,2)
c1.cd(1)
hq_superEarly = TH1F('hq_superEarly', 'hq_all', 400, -10, 20)
hq_range = TH1F( 'hq_range', 'sipm Q, range integral ', 400, -10, 20) # fixed range integral
hq_EarlyRange = TH1F( 'hq_EarlyRange', 'sipm Q, range integral ', 400, -10, 20) # return to baseline integral
hq_flat = TH1F( 'hq_flat', 'sipm Q, range integral ', 400, -10, 20) # return to baseline integral
hq_earlyFlat = TH1F( 'hq_earlyFlat', 'sipm Q, range integral ', 400, -10, 20) # return to baseline integral

hq_lateFlat = TH1F( 'hq_lateFlat', 'sipm Q, range integral ', 400, -10, 20) # return to baseline integral

# getting all data in runlist
for line in runList:

	### CUTS: env/daq level, this is only 55 Volts
	if int(line[jV]) != voltage:
		continue
	if str(line[jSrc]) != source:
		continue

	# getting the path
	path = os.path.join(raw_directory, line[jRun], line[jEv], "PMTtraces.bin")
	data = rb.ReadBlock(path)
	print(path)

	# pulse information
	traces = data.get('traces')
	vScale = data.get('v_scale')
	offset = data.get('v_offset')

	# loop through 100 pulses of one trace object
	# ind loops through pulses
	for ind in range(min(len(traces), maxData)):

		if ind%1000 == 0 and ind > 0:
			c1.Modified()
			c1.Update()
			print("ev#", ind)

		######### Getting Data #########
		# recorded sipm waveform
		s1 = traces[ind, iSipm][0:-30] # shape, last bits are wrong sometimes so skip last 30
		
		s1v = s1 * vScale[ind, iSipm] +  offset[ind, iSipm] # shape x scale + offset (actual pulse data)
		sMin = min(s1v)
		sMax = max(s1v)

		d1 = traces[ind, iDis][0:-30]*vScale[ind, iDis] + offset[ind, iDis] # discriminator

		gap = 1 # pointless to multiply by dt, which is just 1ns
		# gap = data.get('dt')[ind, iSipm] * 10**9
		time = np.array(range(len(s1v))) * gap # + data.get("t0")[ind, iSipm] # ignore the small shift

		######### processing #########
		## ** should really build a processor to handle everything :(

		# integrate over a range (copy zhiheng)
		# apply smoothing filter, try again
		# fit it?

		# handling EARLY pulses:
		#re calculate the charge if max-min is small and/or charge is negative,
		# use the END of the pulse to calculate the waveform 
		# if sMax-sMin < 0.02 and np.argmin(s1v) < 
		#SCRATCH THAT, just skip them for now

		# bline = getBaselineHist(data = s1v, vScale = vScale[ind, iSipm])
		blines = getBaselineAvg(data = s1v, start = 0, end = 200)
		#blinesH = getBaselineHist(data = s1v[0:350], vScale = vScale[ind, iSipm])
		#blinesBack = getBaselineAvg(data = s1v, start = 1250, end = 2000)

		# qRange, rEarly, rLate = minFixRangeIntegrate(data = s1v, bline = blines.get('v'), early = 250, late = 450)
		
		######### CUTTING #########
		earlyFlag = 0 #0 is normal, 1 is early (500ns), 2 is super early (250ns)
		# these pulses are just a flat line starting above 0 with a small area
		if sMax > 0 and sMin > -0.01: #slant cut (recovery pulses, have something above zero with no dip)
			# print("slant")
			continue

		if np.argmin(s1v) < 500 and abs(sMax - sMin) < 0.02: #early pulse cut, which messes up the baseline estimate, so recalculate baseline
			# print("early and short")
			earlyFlag = 1
			blines = getBaselineAvg(data = s1v, start = 0, end = 100) # baseline from the back, only use for small pulses
			if np.argmin(s1v) < 250:
				earlyFlag = 2
			# continue

		##### GET CHARGES #####
		qFull = simpleIntegrate(data = s1v, bline = blines.get('v'))
		# # qFullBack = simpleIntegrate(data = s1v, bline = blinesBack.get('v'))

		qRange, rEarly, rLate = minFixRangeIntegrate(data = s1v, bline = blines.get('v'), early = 250, late = 450)
		# qRangeBack = minFixRangeIntegrate(data = s1v, bline = blinesBack.get('v'), early = 200, late = 300)

		qFlat, fEarly, fLate = toFlatRangeIntegrate(data  = s1v, bline = blines.get('v'), blineRMS = blines.get('e'), extra = 50)

		hq_flat.Fill(qFlat)
		if earlyFlag == 0:
			hq_range.Fill(qRange)
		elif earlyFlag == 1:
			hq_EarlyRange.Fill(qRange)
			hq_earlyFlat.Fill(qFlat)
		elif earlyFlag == 2:
			hq_superEarly.Fill(qFlat)
		if np.argmin(s1v) > 500:
			hq_lateFlat.Fill(qFlat)

		######### Viewing/Troubleshooting #########
		if view == True:
			# for looking at pulses
			ax0.plot(time, s1v, '.', alpha = 0.7)
			ax0.plot(time, d1, '-') # discriminator pulse
			ax0.axvline(rEarly, color = 'k', label = "fixed Range")
			ax0.axvline(rLate, color = 'k')

			ax0.axvline(fEarly, color = 'r', label = "variable Range")
			ax0.axvline(fLate, color = 'r')

			ax0.plot(time, blines.get('v')*np.ones(len(time)), label = "early Baseline")
			# ax0.plot(time, blinesBack.get('v')*np.ones(len(time)))

			print("ev#", ind)
			print( '%0.1E' % blines['v'], '%0.1E' % blines['e'])
			# print('%0.1E' % blinesH[0], '%0.1E' % blinesH[1])
			print(qFull)
			print(qRange)
			print(qFlat)
			print("\n")

			ax0.legend()

			input()
			# plt.pause(0.4)
			ax0.clear()

######### PLOTTING #########
gStyle.SetOptFit(1)
ledge = TLegend(0.6, 0.1, 0.9, 0.4)

# hq_range.SetLineColor(3)
# hq_range.Draw("hist")
# ledge.AddEntry(hq_range, "fixed range", "L")

# hq_EarlyRange.SetLineColor(2)
# hq_EarlyRange.Draw("same")
# ledge.AddEntry(hq_EarlyRange, "early fixed", "L")

hq_flat.SetLineColor(4)
hq_flat.Draw("hist")
ledge.AddEntry(hq_flat, "variable", "L")

hq_earlyFlat.SetLineColor(6)
hq_earlyFlat.Draw("same")
ledge.AddEntry(hq_earlyFlat, "early variable", "L")
ledge.Draw("same")

hq_lateFlat.SetLineColor(3)
hq_lateFlat.Draw("same")
ledge.AddEntry(hq_lateFlat, "late variable", "L")

# hq_superEarly.SetLineColor(1)
# hq_superEarly.Draw("same")
# ledge.AddEntry(hq_superEarly, "super early", "L")
# ledge.Draw("same")

## total hists ##
# c1.cd(2)
# htotal = hq_flat.Clone()
# htotal.Add(hq_earlyFlat)
# htotal.Add(hq_superEarly)
# htotal.Draw("hist")

c1.Modified()
c1.Update()

hFile.cd()
hq_range.Write()
hq_EarlyRange.Write()
hq_flat.Write()
hq_earlyFlat.Write()
hq_superEarly.Write()
# hFile.Close()

# input("Please press enter to continue.")
# final plotting
# ax0.minorticks_on()
# ax0.grid(True)

# ax0.set_ylabel("voltage (mV)")
# ax0.set_xlabel("time (ns)")
# # ax0.set_ylim(-200, 200)

# plt.show()

# exit()




