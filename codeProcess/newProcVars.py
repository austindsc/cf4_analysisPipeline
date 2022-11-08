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
from time import sleep

from sympy import false
sys.path.append('..')
from SBCcode.DataHandling import ReadBinary as rb

from ozTools import *
from ozRCcorrection import *
# from outputManager import *

import ROOT
from ROOT import TFile, TTree
from ROOT import TCanvas, gApplication, gStyle, TLegend, TPaveText
from ROOT import TF1, TH1F, TGraph, TGraphErrors, TFitResult
from ROOT import TSpectrum, TVectorD
from ROOT import TMath, TFormula
from ROOT import gROOT

def on_press(event):
	# print('press', event.key)
	if event.key == 'x' or event.key == 'q':
		exit()
	sys.stdout.flush()


print("Austin's new processor")

SATURATED = 5.
SAMPLING = 16.
BLINEEND = 600.
defRun = 202205271113
defDir = ""
# defDir = "caenTesting"


if __name__ == "__main__":
	pars = argparse.ArgumentParser()
	pars.add_argument("--start", type = int, required = False, default = 0, help = "what event to start in")
	pars.add_argument("--N", type = int, required = False, default = -1, help = "N events to plot/analyze. default -1 runs all")
	pars.add_argument("--run", type = str, required = False, default = defRun, help = "desired run number")
	pars.add_argument("--live", type = bool, required = False, default = False, help = "do viewing with user input")
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
	savePath = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/"+args.outDir
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
	
	if args.live == True or (args.N > 0 and args.N < 1000): # avoid overwriting big files, if user is viewing/testing
		outName = savePath + "LIVE_npv1_" + str(args.run) + ".root"

	hFile = gROOT.FindObject(outName)
	if hFile:
		hFile.Close()
	hFile = TFile(outName, 'RECREATE', 'bprocStaticVars1.py output file')
	print("writing ROOT file to: ", outName, "\n")
	
	### defining TTree and all variables ###
	tt = TTree("tt", "proc static variables, run:" + str(args.run))
	nchans = 3

	tt_evTime = np.empty((1), dtype = "float32")
	tt.Branch("evTime", tt_evTime, "evTime/F")

	tt_blineV = np.empty((nchans), dtype = "float32")
	tt.Branch("blineV", tt_blineV, 'blineV[' + str(nchans) + ']/F')
	tt_blineE = np.empty((nchans), dtype = "float32")
	tt.Branch("blineE", tt_blineV, 'blineE[' + str(nchans) + ']/F')
	tt_min = np.empty((nchans), dtype="float32")
	tt.Branch("minV", tt_min, 'minimum[' + str(nchans) + ']/F')
	tt_xMin = np.empty((nchans), dtype="float32")
	tt.Branch("timeMin", tt_xMin, 'minimum[' + str(nchans) + ']/F')

	tt_qFixed = np.empty((nchans), dtype = "float32")
	tt.Branch("qFixed", tt_qFixed, "qFixed/F")

	tt_qFlat = np.empty((nchans), dtype = "float32")
	tt.Branch("qFlat", tt_qFlat, "qFlat/F")

	tt_tFlatStart = np.empty((nchans), dtype = "float32")
	tt.Branch("tStart", tt_tFlatStart, "tsStart/F")

	tt_tFlatEnd = np.empty((nchans), dtype = "float32")
	tt.Branch("tsEnd", tt_tFlatEnd, "tsEnd/F")

	tt_qRCFixed = np.empty((nchans), dtype = "float32")
	tt.Branch("qRCFixed", tt_qRCFixed, "qRCFixed/F")

	tt_qRCFlat = np.empty((nchans), dtype = "float32")
	tt.Branch("qRCFlat", tt_qRCFlat, "qRCFlat/F")



	### setup plotting ###
	if args.live == True:
		fig, axes = plt.subplots(1, 1, figsize = (8, 6))
		ax0 = axes
		plt.ion()
		plt.show()

	#### LOAD THE DATA	####
	path = os.path.join(raw_directory, str(args.run)+".bin")
	try:
		os.path.exists(path)
	except:
		print("failed to load data")
		exit()
	if not args.caen: #old
		path = os.path.join(raw_directory, run, str(args.g), "PMTtraces.bin")
	
	data = rb.ReadBlock(path)
	traces = data.get("sipm_traces") # [event][channel][waveform value]
	times = data.get("time_stamp")
	# for key, value in data.items(): # inspecting output
	# 	print(key)
	# print(traces.shape)
	# print(traces[0][0])

	#### Event Loop ####
	for ind, [evtime, channels] in enumerate(zip(times[args.start:args.N], traces[args.start:args.N])):
		
		x = np.arange(0, len(channels[0]) * SAMPLING, SAMPLING)
		if(not (channels[0].size == channels[1].size == channels[2].size)):
			print("issue with lengths")
			continue

		if ind%1000 == 0: # printing helper
			print("event", ind, end = "\r")
		if args.live == True:
			print("\nstart index = ", ind+args.start)
		
		#### channel loop
		for ich, wfm in enumerate(channels):

			#### Channel Level Processing ####
			# baseline
			# chargeFixedRange1
			# chargeFixedRange2
			# charge Smart
			# charge RC range 1
			# charge RC range 2
			# charge RC Smart

			# Baselines #
			blines = getBaselineAvg(data = wfm, start = 0, end = int(BLINEEND/SAMPLING))
			tt_blineV[ich] = blines.get("v")
			tt_blineE[ich] = blines.get("e")
			# print("bline Error = ", blines.get("e"))

			# Charge #
			tt_qFixed[ich], rEarly, rLate = minFixRangeIntegrate(data = wfm, bline = blines.get('v'), early = int(180/SAMPLING), late = int(500/SAMPLING))
			tt_qFlat[ich], fEarly, fLate = toFlatRangeIntegrate(data = wfm, bline = blines.get('v'), blineRMS = 5*blines.get('e'), extra = 10)
			tt_tFlatStart[ich] = fEarly*SAMPLING
			tt_tFlatEnd[ich] = fLate*SAMPLING

			# RC Correction
			rcWfm = BaselineRestoration(wfm = wfm-blines.get("v"), tau = 2900, width = SAMPLING)
			tt_xMin[ich] = np.argmin(rcWfm)*SAMPLING
			tt_min[ich] = min(rcWfm)
			tt_qRCFixed[ich],rcEarly, rcLate = minFixRangeIntegrate(data = rcWfm, bline = 0.0, early = int(180/SAMPLING), late = int(500/SAMPLING))
			tt_qRCFlat[ich], rcSEarly, rcSLate = toFlatRangeIntegrate(data  = rcWfm, bline = 0.0, blineRMS = 4*blines.get('e'), extra = 3)
			print("end of RCFlat", rcSLate)
			### plotting channels and printing ###
			colors = ["b", "m", "k"]
			if args.live == True:

				ax0.plot(x, wfm, f'{colors[ich]}', alpha = 0.6, label = f"ch {ich}")
				ax0.axhline(y = blines.get('v'), color = "m", linewidth = 1, alpha = 0.4)

				ax0.axvline(SAMPLING*fEarly, color = f'{colors[ich]}', alpha = 0.6, label = "to baseline limit")
				# ax0.axvline(tt_xMin)
				ax0.axvline(SAMPLING*fLate, color = f'{colors[ich]}', alpha = 0.6)

				ax0.plot(x, rcWfm+blines.get('v'), f'--{colors[ich]}', alpha = 0.8, label = f"rc ch {ich}")
				ax0.axvline(SAMPLING*rcLate, color = f'{colors[ich]}', linestyle = "--", alpha = 0.8, label = "RC limit (fixed)")
				ax0.axvline(SAMPLING*rcEarly, color = f'{colors[ich]}',linestyle = "--", alpha = 0.8)

				ax0.axvline(SAMPLING*(rcSLate), color = f'{colors[ich]}', linestyle = "-.", alpha = 0.8, label = "RC limit (to baseline)")
				ax0.axvline(SAMPLING*rcSEarly, color = f'{colors[ich]}', linestyle = "-.", alpha = 0.8)

				ax0.set_xlim(0, 3000)
				print(f"ch {ich}", end = "\n")
		#### done channel loop

		print("fixed charge", tt_qFixed)
		print("flat charge", tt_qFlat)
		print("tt_qRCFixed", tt_qRCFixed)
		print("tt_qRCFlat", tt_qRCFlat)

		print("tt_min", tt_min)
		
		#### Event level Processing ####
		tt_evTime[0] = evtime

		#### Fill Tree ####
		tt.Fill()

		##### plotting event level #####
		if args.live == True:

			# print("event time %i" % evtime)

			plt.pause(0.05)
			ax0.minorticks_on()
			ax0.grid(True)
			ax0.legend()
			plt.draw()
			plt.pause(0.1)
			
			while(plt.waitforbuttonpress(timeout = -1) == False):
				plt.pause(0.1)
			plt.connect('key_press_event', on_press)
			ax0.clear()


	#### Closing and Writing #####
	hFile.cd()
	tt.Write()
	hFile.Close()

	print("Wrote output to", outName)





