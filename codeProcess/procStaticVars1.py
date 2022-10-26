#!/usr/local/bin/python3

"""
	version1 of code to process the data, and output the variables to a TTree for easy analysis later

	fitting is currently commented out and not written to TTree (save time)

running instructions

	python3 procStaticVars1.py --run 20210729_4
		this will fully process all the events in given run and output "outProcStatic1_20210729_4.root" to savePath variable

	python3 procStaticVars.py --help
		for other options explained
"""

# need to add lines where all variables are calculated from
# maybe difference betweeen data is becaus of baselines

import numpy as np
import os
import matplotlib.pyplot as plt
import argparse
from array import array

import sys
sys.path.append('..')
from SBCcode.DataHandling import ReadBinary as rb

from ozTools import *

import ROOT
from ROOT import TFile, TTree
from ROOT import TCanvas, gApplication, gStyle, TLegend, TPaveText
from ROOT import TF1, TH1F, TGraph, TGraphErrors, TFitResult
from ROOT import TSpectrum, TVectorD
from ROOT import TMath, TFormula
from ROOT import gROOT

# must change filename when using 'live'

print("Austin's makeshift static variable processor")

# index of recorded waveforms in the event
iSipm = 0 # sipm data is first
iDis = 1 # discriminator data second

###### parsing ######
pars = argparse.ArgumentParser()
pars.add_argument("--first", type = int, required = False, default = 0, help = "what event to start at")
pars.add_argument("--g", type = int, required = False, default = 0, help = "what pulse group to start in (folder)")
pars.add_argument("--N", type = int, required = False, default = -1, help = "N events to loop through. Default -1 means process all")
pars.add_argument("--run", type = str, required = False, default = "x_20210726_0", help = "randomly selected default run")
pars.add_argument("--live", type = int, required = False, default = False, help = "view each pulse, with user pressing Enter")
pars.add_argument("--view", type = int, required = False, default = False, help = "view final plots")
pars.add_argument("--outName", type = str, required = False, help = "if custom naming of output file is desired, input name here (without .root suffix)")
args=pars.parse_args()
run = args.run
live = bool(args.live)
view = bool(args.view)

##### output and input file #####
savePath = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/"
raw_directory = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/rawdata/"
outName = savePath + "outProcStatic1_" + str(run) + ".root"
if args.outName != None:
	outName = str(savePath + args.outName + str(run) + ".root")

if args.live == True or args.view == True or (args.N > 0 and args.N < 1000): # avoid overwriting big files, if user is viewing/testing
	outName = savePath + "VIEW_outProcStatic1_" + str(run) + ".root"

hFile = gROOT.FindObject(outName)
if hFile:
	hFile.Close()
hFile = TFile(outName, 'RECREATE', 'bprocStaticVars1.py output file')
print("writing ROOT file to: ", outName, "\n")

### defining TTree and all variables ###
tt = TTree("tt", "proc static variables, run:" + str(run))

blinePre = 400
# 500 for mega trigger
# 40 for medium gas trigger

height1 = np.empty((1), dtype = "float32")
tt.Branch("hbline_min", height1, "hbline_min/F")

blineVal = np.empty((1), dtype = "float32")
tt.Branch("blineV", blineVal, "blineV/F")

blineE = np.empty((1), dtype = "float32")
tt.Branch("blineE", blineE, "blineE/F")

qr1 = np.empty((1), dtype = "float32")
tt.Branch("qr300", qr1, "qr300/F")

qr2 = np.empty((1), dtype = "float32")
tt.Branch("qr200", qr2, "qr200/F")

qr3 = np.empty((1), dtype = "float32")
tt.Branch("qr800", qr3, "qr800/F")

qs1 = np.empty((1), dtype = "float32")
tt.Branch("qs50", qs1, "qs50/F")

tstart = np.empty((1), dtype = "float32")
tt.Branch("tStart", tstart, "tsStart/F")

tend = np.empty((1), dtype = "float32")
tt.Branch("tsEnd", tend, "tsEnd/F")

# f1q = np.empty((1), dtype = "float32")
# tt.Branch("f1q", f1q, "f1q/F")

# f1r = np.empty((1), dtype = "float32")
# tt.Branch("f1r", f1r, "f1r/F")

# f1f = np.empty((1), dtype = "float32")
# tt.Branch("f1f", f1f, "f1f/F")

# f1chi2 = np.empty((1), dtype = "float32")
# tt.Branch("f1chi2", f1chi2, "f1chi2/F")


######### BUILDING charge HISTS #########
qMax = 20
hq_fAll = TH1F("hq_fAll", "hq_fAll", 200, -1, qMax)

hq_fChi2 = TH1F("hq_fChi2", "hq_fChi2", 200, -1, qMax)

hq_range1 = TH1F("hq_range1", "hq_range1", 200, -1, qMax)
hq_range2 = TH1F("hq_range2", "hq_range2", 200, -1, qMax)
hq_range3 = TH1F("hq_range3", "hq_range3", 200, -1, qMax)
hq_smart1 = TH1F("hq_smart1", "hq_smart1", 200, -1, qMax)

hTest = TH1F("hTest", "hTest", 200, -1, qMax)

hh_static = TH1F("hh_static", "hh_static", 200, 0, 0.1)
hh_fits = TH1F("hh_fits", "hh_fits", 200, 0, 0.1)

fitOpt = "WQN0"
specOpt = "goff"
if live == True: #setup plotting
	### setup plotting ###
	fig, axes = plt.subplots(1, 1, figsize = (8, 6))
	ax0 = axes
	plt.ion()

	# c1 = TCanvas('c1', 'c1', 600, 600)
	fitOpt = "WN0"
	specOpt = "goff"


##### pulse group loop #####
totalN = 0
group = args.g
path = os.path.join(raw_directory, run, str(group), "PMTtraces.bin")
print(path)
while os.path.exists(path) == True:

	## loading data
	data = rb.ReadBlock(path)
	group = group + 1
	path = os.path.join(raw_directory, run, str(group), "PMTtraces.bin")

	# data is a dictionairy, you can print all its members (google python dictionairy)
	# getting pulse information
	traces = data.get('traces')
	vScale = data.get('v_scale')
	offset = data.get('v_offset')
	
	# loop through pulses (aka events), loops through all events if args.N == 0
	# for ind in range(args.first, (lambda: (args.first + args.N, len(traces))[args.N > len(traces) or args.N == 0])()):
	for ind in range(args.first, len(traces)):
		totalN = totalN + 1
		if totalN > args.N and args.N != -1:
			break
		if ind%1000 == 0: # printing helper
			print("at event", ind)

		# recorded sipm waveform
		s1 = traces[ind, iSipm][0:-100] # shape, last few values are corrupt or something
		s1v = s1 * vScale[ind, iSipm] +  offset[ind, iSipm] # shape x scale + offset
		d1 = traces[ind, iDis] # discriminator
		# # trigger from discriminator
		# trigger = traces['traces'][ind,1] * traces['v_scale'][ind,1] + traces['v_offset'][ind,1]

		gap = 1 # pointless to multiply by dt, which is just 1ns
		# gap = data.get('dt')[ind, iSipm] * 10**9
		time = np.array(range(len(s1v))) * gap #+ data.get("t0")[ind, iSipm]

		######## doing processing ########
		minTime = np.argmin(s1v)
		blines = getBaselineAvg(data = s1v, start = 0, end = blinePre)
		earlySpread = blines.get('v') - min(s1v[0:40]) # may want to cut on this variable later (at spread < around 0.007)

		if(earlySpread > 0.01):
			continue

		# if minTime < 300: # pulse is too early
		# 	blines = getBaselineAvg(data = s1v, start = minTime+600, end = minTime+900)
		# 	if live == True:
		# 		print("using end as baseline")

		qRange, rEarly, rLate = minFixRangeIntegrate(data = s1v, bline = blines.get('v'), early = 150, late = 300)
		qRange2, rEarly2, rLate2 = minFixRangeIntegrate(data = s1v, bline = blines.get('v'), early = 100, late = 200)
		qRange3, rEarly3, rLate3 = minFixRangeIntegrate(data = s1v, bline = blines.get('v'), early = 200, late = 900)
		qSmart1, fEarly, fLate = toFlatRangeIntegrate(data  = s1v, bline = blines.get('v'), blineRMS = blines.get('e'), extra = 10)

		# building TH1 (and flipping)
		hwf = TH1F("hwf", "hwf", len(s1v), time[0], time[-1])
		for i in range(len(s1v)):
			hwf.SetBinContent(i+1, -s1v[i])

		# if live == True:
		# 	hwf.Draw("hist") # root drawing first so spectrum peaks show up

		# peak finding (currently unused)
		spec = TSpectrum(1)
		spec.SetAverageWindow(6)
		ROOT.gVerboseLevel = ROOT.kError
		ROOT.gErrorIgnoreLevel = ROOT.kFatal
		spec.Search(hwf, 10., specOpt + "noMarkov", 0.5)

		###### fitting functions ######
		flat1 = TF1("flat1", "pol0", 0, 250)
		flat1.SetParameters(-blines.get('v'))

		theight1 = blines.get('v') - min(s1v)
		height2 = spec.GetPositionY()[0] + blines.get('v')

		hh_static.Fill(theight1)
		hh_fits.Fill(height2)

		# 0 - area
		# 1 - rise time
		# 2 - fall time (Tl)
		# 3 - shift
		# 4 - baseline
		"""
		pulse1 = TF1("pulse1", "[0]/[2] * (TMath::Exp(-(x-[3])/([2] + [1])) - TMath::Exp(-(x - [3])/[1])) + [4]", spec.GetPositionX()[0] - 300, spec.GetPositionX()[0] + 800)
		pulse1.SetParameters(spec.GetPositionY()[0]*200, 70, 80, minTime, -blines.get('v'))
		pulse1.SetParNames("area", "riseT", "fallT", "shift", "endBline")
		pulse1.SetParLimits(0, 0.5, 100)
		pulse1.SetParLimits(1, 10, 100)
		pulse1.SetParLimits(2, 40, 250)
		pulse1.SetParLimits(3, minTime - 150, minTime + 100)
		# pulse1.SetParLimits(4, -blines.get('v'), -blines.get('v'))
		pulse1.SetParLimits(4, -blines.get('v'), -blines.get('v'))

		hwf.Fit(pulse1, fitOpt, "goff", spec.GetPositionX()[0] - 100, spec.GetPositionX()[0] + 600)
		"""

		####### fill TTreee #######
		blineVal[0] = float(blines.get('v'))
		height1[0] = blines.get('v') - min(s1v)
		qr1[0] = qRange
		qr2[0] = qRange2
		qr3[0] = qRange3
		qs1[0] = qSmart1
		tend[0] = fLate
		tstart[0] = fEarly
		
		"""
		f1q[0] = pulse1.GetParameter(0)
		f1r[0] = pulse1.GetParameter(1)
		f1f[0] = pulse1.GetParameter(2)
		f1chi2[0] = pulse1.GetChisquare()
		"""

		tt.Fill() #adding a new row to the ttree (saving the data)

		###### filling histograms, TTree with values #####
		# qFit = pulse1.GetParameter(0)
		# chi2Fit = pulse1.GetChisquare()

		# hq_fAll.Fill(qFit)
		# if chi2Fit <= 2.0:
		# 	hq_fChi2.Fill(qFit)

		## filling simple histograms
		hq_range1.Fill(qRange)
		hq_range2.Fill(qRange2)
		hq_range3.Fill(qRange3)
		hq_smart1.Fill(qSmart1)

		if live == True:
			# root drawing
			# hwf.Draw("hist")
			# pulse1.Draw("same")
			# flat1.Draw("hist same")

			# c1.Modified()
			# c1.Update()

			#hold for user input
			print("\n event #", ind)
			print("minTime = ", minTime)
			print("peaks found = ", spec.GetPositionX()[0])
			print("blines = ", blines)
			print("qRange = ", qRange)
			print("qSmart = ", qSmart1)
			# print("Qfit = ", qFit, "    ch2 = ", chi2Fit, "   NDF = ", pulse1.GetNDF())
			# del pulse1

			ax0.clear() # somehow reset but keep last N pulses
			plt.pause(0.1)
			ax0.hlines(blines.get('v'), time[0], time[-1])
			ax0.plot(time, s1v, '.', alpha = 0.5)
			print(s1v[-10:-1])
			input()



		del hwf

hFile.cd()
tt.Write()

num = TVectorD(1)
num[0] = blinePre
num.Write("blinePre")
# to access do:
# TVectorD * a = (TVectorD*) _file0->Get("blinePre")
# (*a)[0]


if live == False and view == True:
	c1 = TCanvas('c1', 'c1', 600, 400)

	hq_smart1.SetLineColor(1) # black
	hq_smart1.Draw("hist")

	hq_range1.SetLineColor(2) # red
	hq_range1.Draw("hist same")

	hq_range2.SetLineColor(3) # green
	hq_range2.Draw("hist same")

	hq_range3.SetLineColor(4) # blue
	hq_range3.Draw("hist same")

	# hq_fAll.SetLineColor(4) # blue
	# hq_fAll.Draw("hist same")

	# hq_fChi2.SetLineColor(6) # pink
	# hq_fChi2.Draw("hist same")

	c1.Update()

	# c2 = TCanvas("c2", "c2", 600, 400)

	# hh_static.Draw("hist")
	# hh_fits.SetLineColor(3)
	# hh_fits.Draw("hist same")
	# c2.Update()

	ROOT.gApplication.Run()

hFile.Close()

