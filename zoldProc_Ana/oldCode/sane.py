#!/usr/local/bin/python3

#stock python stuff
import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
import time as chrono

## this is pyroot
from ROOT import ROOT
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
# 1. load in all runs (date and number)
# 2. get N events? or just try until read fails
# 3. search 

# can I remove noise from pulses by looking at noise in the disc. signal??

#### plans? ####
# okay lets look at SPE distributions for small runs and big runs

#compare low PE spectra to get better QPE values

# goal: find what 1PE is
# maybe we need to apply a moving average/filter

# also try looking at rates of low Q pulses (if rate increases with voltage, not noise)
#########################

######### PARSING #########
pars = argparse.ArgumentParser()
pars.add_argument("--view", type = int, required = False, default = False)
pars.add_argument("--max", type = int, required = False, default = 100000)
pars.add_argument("--skip", type = int, required = False, default = 1)
args = pars.parse_args()

maxData = args.max
skip = args.skip
view = bool(args.view) # true (1) for viewing pulses

######### RUN LISTS, DATA INDICES, DIRECTORIES #########
runPath = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/atHome/runLists/"
runFile = "liquidBig.txt"

voltage = 55
source = 'Cs137'
ampGain = 20
ddate = "20201116"

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

# # root files
hFile = gROOT.FindObject('SPE1.root')
if hFile:
	hFile.Close()
hFile = TFile('SPE1.root', 'RECREATE', 'speAna1 test file')


######### BUILDING HISTS #########
h2_range = TH2F('h2_range', 'h2_range', 400, 0, 1500, 100, -2, 10)
h2_flat = TH2F('h2_flat', 'h2_range', 400, 0, 1500, 100, -2, 10)
h2_height = TH2F("h2_height", "h2_height", 400, 0, 1500, 1000, -0.01, 0.3)

hR_all = TH1F('hR_all', 'hq_all', 200, -10, 20)
hR_noSlant = TH1F('hR_noSlant', 'sipm Q, range integral ', 200, -10, 20) # fixed range integral

hF_all = TH1F('hF_all', 'hq_all', 2000, -10, 150)
hF_noSlant = TH1F('hF_noSlant', 'sipm Q, range integral ', 2000, -10, 150) # fixed range integral

hH_all = TH1F('hH_all', 'hH_all', 1000, -0.01, 0.3)

######### looping on data #########=
for line in runList: # getting all data in runlist

	### CUTS: env/daq level, voltage, source, etc
	if int(line[jV]) != voltage:
		continue
	if str(line[jSrc]) != source:
		continue
	if str(line[jRun][:8]) != ddate:
		continue
	if int(line[jGain]) != ampGain:
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
	print("N events =", len(traces))
	for ind in range(0, min(len(traces), maxData), skip):

		if ind%1000 == 0 and ind > 0:
			print("ev#", ind)

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

		######### processing #########

		# simple blines calc (FIRST 100 pts)
		blines = getBaselineAvg(data = s1v, start = 0, end = 200)
		qRange, rEarly, rLate = minFixRangeIntegrate(data = s1v, bline = blines.get('v'), early = 250, late = 450)
		qFlat, fEarly, fLate = toFlatRangeIntegrate(data  = s1v, bline = blines.get('v'), blineRMS = blines.get('e'), extra = 50)
		qFlat2, fEarly2, fLate2 = toFlatRangeIntegrate(data  = s1v, bline = blines.get('v'), blineRMS = blines.get('e'), extra = 250)

		minTime = np.argmin(s1v)
		height = minAverage(data = s1v, bline = blines.get('v'), spread = 8)
		
		######### CUTTING and FILLING #########		
		hR_all.Fill(qRange)
		hF_all.Fill(qFlat)
		hH_all.Fill(height)

		h2_range.Fill(minTime, qRange)
		h2_flat.Fill(minTime, qFlat)
		h2_height.Fill(minTime, height)

		if np.argmin(s1v) < 300:
			if view == True:
				print("too early \n")
			continue

		# these pulses are just a flat line starting above 0 with a small area
		if sMax > 0 and sMin > -0.01 or (sMax - sMin) < 0.008: #slant cut (recovery pulses, have something above zero with no dip)
			if view == True:
				print("SLANT \n")
			continue

		hR_noSlant.Fill(qRange)
		hF_noSlant.Fill(qFlat)

		######### Viewing/Troubleshooting #########
		if view == True:
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
			print(blines.get('v'), sMax - sMin)
			# print( '%0.1E' % blines['v'], '%0.1E' % blines['e'])
			# print('%0.1E' % blinesH[0], '%0.1E' % blinesH[1])
			print(qRange)
			print(qFlat)
			print(qFlat2, "\n")

			input()
			# plt.pause(0.4)
			ax0.clear()

######### PLOTTING #########
# c1 = TCanvas("c1", "c1", 600, 600)
# c1.Divide(1,2)
# c1.cd(1)

# gStyle.SetOptFit(1)
# ledge = TLegend(0.6, 0.1, 0.9, 0.4)

# # all, range
# hR_all.SetLineColor(1)
# hR_all.Draw("hist")
# ledge.AddEntry(hR_all, "all - fixed", "L")

# #all, flat
# hF_all.SetLineColor(2)
# hF_all.Draw("same")
# ledge.AddEntry(hF_all, "all - variable", "L")

# #no slant, range
# hR_noSlant.SetLineColor(3)
# hR_noSlant.Draw("same")
# ledge.AddEntry(hR_noSlant, "no slants - fixed", "L")

# # no slant, flat
# hF_noSlant.SetLineColor(4)
# hF_noSlant.Draw("same")
# ledge.AddEntry(hF_noSlant, "no slants - variable", "L")

# ledge.Draw("same")

# c1.cd(2)
# ledge2 = TLegend(0.6, 0.1, 0.9, 0.4)

# hH_all.Draw("hist")
# ledge2.AddEntry(hH_all, "pulse height", "L")

# c1.Modified()
# c1.Update()


# c2 = TCanvas("c2", "c2", 600, 600)
# c2.Divide(1,3)

# c2.cd(1)
# h2_range.Draw("colz")

# c2.cd(2)
# h2_flat.Draw("colz")

# c2.cd(3)
# h2_height.Draw("colz")

# c2.Modified()
# c2.Update()

hFile.cd()

hR_all.Write()
hR_noSlant.Write()

hF_all.Write()
hF_noSlant.Write()

hH_all.Write()

hFile.Close()

# input("Please press enter to continue.")
# final plotting
# ax0.minorticks_on()
# ax0.grid(True)

# ax0.set_ylabel("voltage (mV)")
# ax0.set_xlabel("time (ns)")
# # ax0.set_ylim(-200, 200)

# plt.show()

# exit()




