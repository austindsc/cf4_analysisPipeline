#!/usr/local/bin/python3.9

"""
	use this for viewing waveforms. Minimal processing done here. good 'skeleton' for viewing data
"""

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

print("Austin's new processor")

defRun = 202205271113
defDir = ""

if __name__ == "__main__":
	pars = argparse.ArgumentParser()
	pars.add_argument("--first", type = int, required = False, default = 0, help = "what event to start in")
	pars.add_argument("--N", type = int, required = False, default = 10, help = "N events to plot/analyze")
	pars.add_argument("--run", type = str, required = False, default = defRun, help = "desired run number")
	pars.add_argument("--live", type = bool, required = False, default = False, help = "do viewing with user input")
	pars.add_argument("--view", type = int, required = False, default = False, help = "view final plots")
	pars.add_argument("--caen", type = bool, required = False, default = True, help = "is this new data (on caen)")
	pars.add_argument("--outDir", type = str, required = False, default = defDir, help = "name of ouptut directory in main directory")
	pars.add_argument("--outName", type = str, required = False, help = "if custom naming of output file is desired, input name here (without .root suffix)")
	
	## for old data type (non CAEN)
	pars.add_argument("--g", type = int, required = False, default = 0, help = "what pulse group to start on. not for CAEN")
	# index of recorded waveforms in the event
	iSipm = 0 # sipm data is first
	iDis = 1 # discriminator data second
	args=pars.parse_args()
	
	# files and paths
	raw_directory = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/rawData/caen"
	savePath = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/"+args.outDir+"/"
	if not os.path.exists(savePath):
		print("savePath doesn't exist")
		print(savePath)
		exit()
	fname = str(args.run)
	
	if not args.caen:
		raw_directory = "old path"
		fname = "other default file and path"
	
	outName = str(savePath + "npv1_"+ fname +".root")
	if args.outName != None:
		outName = str(savePath + "npv1_"+ args.outName + "_" + fname + ".root")
	
	if args.live == True or args.view == True or (args.N > 0 and args.N < 1000): # avoid overwriting big files, if user is viewing/testing
		outName = savePath + "VIEW_npv1_" + str(args.run) + ".root"

	hFile = gROOT.FindObject(outName)
	if hFile:
		hFile.Close()
	hFile = TFile(outName, 'RECREATE', 'bprocStaticVars1.py output file')
	print("writing ROOT file to: ", outName, "\n")
	
	### defining TTree and all variables ###
	tt = TTree("tt", "proc static variables, run:" + str(args.run))
	nchans = 3
	height1 = np.empty((nchans), dtype = "float32")
	tt.Branch("hbline_min", height1, "hbline_min/F")

	blineVal = np.empty((nchans), dtype = "float32")
	tt.Branch("blineV", blineVal, "blineV/F")

	blineE = np.empty((nchans), dtype = "float32")
	tt.Branch("blineE", blineE, "blineE/F")

	qr1 = np.empty((nchans), dtype = "float32")
	tt.Branch("qr300", qr1, "qr300/F")

	qr2 = np.empty((nchans), dtype = "float32")
	tt.Branch("qr200", qr2, "qr200/F")

	qr3 = np.empty((nchans), dtype = "float32")
	tt.Branch("qr800", qr3, "qr800/F")

	qs1 = np.empty((nchans), dtype = "float32")
	tt.Branch("qs50", qs1, "qs50/F")

	tstart = np.empty((nchans), dtype = "float32")
	tt.Branch("tStart", tstart, "tsStart/F")

	tend = np.empty((nchans), dtype = "float32")
	tt.Branch("tsEnd", tend, "tsEnd/F")

	### setup plotting ###
	if args.live == True:
		fig, axes = plt.subplots(1, 1, figsize = (8, 6))
		ax0 = axes
		plt.ion()
		plt.show()
		plt.plot()
		

	## LOAD THE DATA	
	path = os.path.join(raw_directory, str(args.run)+".bin")
	try:
		os.path.exists(path+"zz")
	except:
		print("failed")
	exit()
	if not args.caen: #old
		path = os.path.join(raw_directory, run, str(args.g), "PMTtraces.bin")
	
	data = rb.ReadBlock(path)
	for key, value in data.items():
		print(key)
	# print(data.get("sample_rate"))
	
	traces = data.get("sipm_traces") # [event][channel][waveform value]
	times = data.get("time_stamp")
	# print(traces.shape)
	# print(traces[0][0])

	# loop through 1000 pulses (or less if not many events)
	for ind, [time, channels] in enumerate(zip(times[0:args.N], traces[0:args.N])):
		
		print("event time %i" % time)
		
		# channel loop
		for ich, wfm in enumerate(channels):
			x = np.arange(0, len(wfm) * 16, 16)
			
			if args.live:
				ax0.plot(x, wfm, alpha = 0.4, label = "ch 0")
			


		if ind%100 == 0: # printing helper
			print("event", ind)



		######## do processing ########
		# minTime = np.argmin(s1v)
		# blines = getBaselineAvg(data = s1v, start = 0, end = 200)
		# 
		# if minTime < 300: # pulse is too early
		# 	blines = getBaselineAvg(data = s1v, start = minTime+400, end = minTime+900)
		# 	if args.live == True:
		# 		print("using end as baseline")
		# 
		# qRange, rEarly, rLate = minFixRangeIntegrate(data = s1v, bline = blines.get('v'), early = 250, late = 450)
		# qSmart1, fEarly, fLate = toFlatRangeIntegrate(data  = s1v, bline = blines.get('v'), blineRMS = blines.get('e'), extra = 0)
		# 

	# 	if args.live == True:
	# 		ax0.clear() # somehow reset but keep last N pulses
	# 		ax0.plot(x0, tr0, alpha = 0.4, label = "ch 0")
	# 		ax0.plot(x1, tr1, alpha = 0.4, label = "ch 1")
	# 		ax0.plot(x2, tr2, alpha = 0.4, label = "ch 2")
	# 		#hold for user input
	# 		# print("\n event #", ind)
	# 		# print("minTime = ", minTime)
	# 		# print("blines = ", blines)
	# 		# print("qRange = ", qRange)
	# 		# print("qSmart = ", qSmart1)
	# 		plt.pause(0.1)
	# 		ax0.minorticks_on()
	# 		ax0.grid(True)
	# 		ax0.legend()
	# 		plt.show()
	# 		input("\n press a key to see next wfm")
	# 
	# if not args.live:
	# 	ax0.minorticks_on()
	# 	ax0.grid(True)
	# 	ax0.legend()
	# 	# ax0.set_ylim(-200, 200)
	# 
	# 	plt.show()





