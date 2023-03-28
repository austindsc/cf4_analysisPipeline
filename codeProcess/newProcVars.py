
import numpy as np
import os
import matplotlib.pyplot as plt
import argparse
from array import array

import sys
from time import sleep

# from sympy import false
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

## must run with python 3.9 for multiprocessing to work ##

def on_press(event):
	# print('press', event.key)
	if event.key == 'x' or event.key == 'q' or event.key == 'escape':
		exit()
	sys.stdout.flush()


	#	! hello
	# TODO yoo
	
"""
1. measure how long the pulse is saturated
2. get a better 'width' estimator 
3. don't skip everything. just fix the error from the bad average in findFlat

4. then we can use 1 and 2 to 'correct' for saturation by using charge vs height and width, 
and the RC correction?
^^do this in postprocessing


other
1. add Nmsising and other stuff to print output
2. don't do any cutting here. do it in analysis (even though i want to here!) ?

try adding a fitter? No. too slow. not what we want. time to move forward.

"""

"""
with ken Nov 11

fixed/toFlat integral --> RC correction integral
FWHM on 'good pulses' --> RC charge or fixed charge
maybe 'to flat' estimator is better - doesn't need pulse height

overshoot --> linear with pulse size

Nsaturation - number of bins at minimum

"""


"""
__ what does this code do? __
calculate fixed charge from minimum
no pulse (small fixed charge or min-baseline is small)
	- skip
early pulse
	- do normal (remove early pulses in postprocessing)
saturated
	- count number of bins at maximum (going left to right)
normal pulse
	- fixed, toFlat charge
	- RC correction, fixed and toFlat charge


"""

print("Austin's new processor")

SATURATED = 6.
SAMPLING = 16.
BLINEEND = 500.
NULLVAL = 0.0
SKIPENDTIME = 500.
defRun = 202205271113
defDir = "2022/05"
# defDir = "caenTesting"


if __name__ == "__main__":
	pars = argparse.ArgumentParser()
	pars.add_argument("--start", type = int, required = False, default = 0, help = "what event to start in")
	pars.add_argument("--N", type = int, required = False, default = -1, help = "N events to plot/analyze. default -1 runs all")
	pars.add_argument("--run", type = str, required = False, default = defRun, help = "desired run number")
	pars.add_argument("--inDir", type = str, required = False, default = defDir, help = "path from rawData to look for run")

	pars.add_argument("--live", type = bool, required = False, default = False, action=argparse.BooleanOptionalAction, help = "do viewing with user input")
	pars.add_argument("--caen", type = bool, required = False, default = True, help = "is this new data (on caen)")
	pars.add_argument("--outDir", type = str, required = False, default = "", help = "name of ouptut directory in main directory")
	pars.add_argument("--outName", type = str, required = False, help = "if custom naming of output file is desired, input name here (without .root suffix)")
	
	## for old data type (non CAEN)
	pars.add_argument("--g", type = int, required = False, default = 0, help = "what pulse group to start on. not for CAEN")
	# index of recorded waveforms in the event
	iSipm = 0 # sipm data is first
	iDis = 1 # discriminator data second
	args=pars.parse_args()
	
	# files and paths
	raw_directory = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/rawData/"+args.inDir
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
	
	if args.live == True or (args.N > 0 and args.N < 900): # avoid overwriting big files, if user is viewing/testing
		outName = savePath + "LIVE_npv1_" + str(args.run) + ".root"


	### Used to open root file here
	
	### defining TTree and all variables ###
	tt = TTree("tt", "proc static variables, run:" + str(args.run))
	nchans = 3

	tt_chInfo = np.zeros((nchans), dtype = "int")
	tt.Branch("chInfo", tt_chInfo, "chInfo[" + str(nchans) +"]/I")
	pulseDict = {"early":-1, "empty": 0, "normal": 1, "saturated":2	}

	tt_evTime = np.zeros((1), dtype = "float32")
	tt.Branch("evTime", tt_evTime, "evTime/F")

	tt_blineV = np.zeros((nchans), dtype = "float32")
	tt.Branch("blineV", tt_blineV, 'blineV[' + str(nchans) + ']/F')
	tt_blineE = np.zeros((nchans), dtype = "float32")
	tt.Branch("blineE", tt_blineV, 'blineE[' + str(nchans) + ']/F')
	tt_min = np.zeros((nchans), dtype="float32")
	tt.Branch("minV", tt_min, 'minV[' + str(nchans) + ']/F')
	##add 3 bin averag around minimum
	tt_xMin = np.zeros((nchans), dtype="float32")
	tt.Branch("timeMin", tt_xMin, 'timeMin[' + str(nchans) + ']/F')

	tt_qFixed = np.zeros((nchans), dtype = "float32")
	tt.Branch("qFixed", tt_qFixed, "qFixed[" + str(nchans) + "]/F")
	tt_qFlat = np.zeros((nchans), dtype = "float32")
	tt.Branch("qFlat", tt_qFlat, "qFlat[" + str(nchans) + "]/F")
	tt_tFlatStart = np.zeros((nchans), dtype = "float32")
	tt.Branch("tFlatStart", tt_tFlatStart, "tFlatStart[" + str(nchans) + "]/F")
	tt_tFlatEnd = np.zeros((nchans), dtype = "float32")
	tt.Branch("tFlatEnd", tt_tFlatEnd, "tFlatEnd[" + str(nchans) + "]/F")
	## write better FWHM width estimator
	## Nminimum bins 
	tt_qRCFixed = np.zeros((nchans), dtype = "float32")
	tt.Branch("qRCFixed", tt_qRCFixed, "qRCFixed[" + str(nchans) + "]/F")
	tt_qRCFlat = np.zeros((nchans), dtype = "float32")
	tt.Branch("qRCFlat", tt_qRCFlat, "qRCFlat[" + str(nchans) + "]/F")

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
	# print("maxTime = ", max(times))
	# for key, value in data.items(): # inspecting output
	# 	print(key, value)
	# print(traces.shape)
	# print(traces[0][0])

	NSaturated = 0
	Nmissing = 0
	Nearly = 0

	#### Event Loop ####
	Nevents = len(times[args.start:args.N])
	for ind, [evtime, channels] in enumerate(zip(times[args.start:args.N], traces[args.start:args.N])):
		
		x = np.arange(0, len(channels[0]) * SAMPLING, SAMPLING)
		if(not (channels[0].size == channels[1].size == channels[2].size)):
			print("issue with lengths")
			continue

		pscale = 1000
		if ind%pscale == 1: # printing helper
			print("event %ik/%ik" %(ind/pscale, Nevents/pscale), end = "\r")
		if args.live == True:
			print("\nstart index = ", ind+args.start)
		
		#### channel loop
		for ich, wfm in enumerate(channels):

			#### Channel Level Processing ####
			# Baselines #
			# blines = getBaselineAvg(data = wfm, start = 0, end = int(BLINEEND/SAMPLING))
			blines = getBaselineMode(data = wfm, start = 0, end = int(BLINEEND/SAMPLING))
			tt_blineV[ich] = blines.get("v")
			tt_blineE[ich] = blines.get("e")

			binMin = np.argmin(wfm[0:-int(SKIPENDTIME/SAMPLING)])
			vMin = min(wfm[0:-int(SKIPENDTIME/SAMPLING)])
			tt_xMin[ich] = binMin*SAMPLING
			tt_min[ich] = vMin

			# Charge #
			testqFlat, rEarly, rLate = minFixRangeIntegrate(data = wfm, bline = blines.get('v'), early = int(180/SAMPLING), late = int(500/SAMPLING), iMin = binMin)
			tt_qFixed[ich] = testqFlat
			tt_chInfo[ich] = pulseDict.get("normal")
			fEarly, fLate = 0, 0
			rcEarly, rcLate = 0, 0
			rcSEarly, rcSLate = 0, 0
			tt_qFlat[ich] = NULLVAL
			tt_tFlatStart[ich] = NULLVAL
			tt_tFlatEnd[ich] = NULLVAL
			tt_qRCFixed[ich] = NULLVAL
			tt_qRCFlat[ich] = NULLVAL

			height = abs(vMin - blines.get("v"))
			
			if height < 5: # zero pulse, or too small
				# print(f"ch{ich} EMPTY")
				Nmissing+=1
				tt_chInfo[ich] = pulseDict.get("empty")
				pass
			elif binMin < int(BLINEEND/SAMPLING):
				# print(f"ch{ich} EARLY")
				# print("binMin x sample = ", binMin * SAMPLING)
				Nearly +=1
				tt_chInfo[ich] = pulseDict.get("early")
			elif vMin <= SATURATED:
				# print(f"ch{ich} SATURATED")
				tt_chInfo[ich] = pulseDict.get("saturated")
				NSaturated += 1
				## do saturated processing ##
			else:
				tt_qFlat[ich], fEarly, fLate = toFlatRangeIntegrate(data = wfm, bline = blines.get('v'), blineRMS = 4*blines.get('e'), extra = 5, iMin = binMin)
				tt_tFlatStart[ich] = fEarly*SAMPLING
				tt_tFlatEnd[ich] = fLate*SAMPLING

				# RC Correction
				rcWfm = BaselineRestoration(wfm = wfm-blines.get("v"), tau = 2900, width = SAMPLING)
				tt_qRCFixed[ich],rcEarly, rcLate = minFixRangeIntegrate(data = rcWfm, bline = 0.0, early = int(180/SAMPLING), late = int(700/SAMPLING), iMin = binMin)
				tt_qRCFlat[ich], rcSEarly, rcSLate = toFlatRangeIntegrate(data  = rcWfm, bline = 0.0, blineRMS = 4*blines.get('e'), extra = 5, iMin = binMin)

			### plotting channels and printing ###
			colors = ["b", "m", "k"]
			if args.live == True:
				
				ax0.plot(x, wfm, f'{colors[ich]}', alpha = 0.6, label = f"ch {ich}")
				ax0.axhline(y = blines.get('v'), color = "m", linewidth = 1, alpha = 0.4)

				if tt_chInfo[ich] == 1:
					ax0.axvline(SAMPLING*fEarly, ymax = 1., ymin = 0.5, color = f'{colors[ich]}', linestyle = '-', label = "normal - flat")
					# ax0.axvline(tt_xMin)
					ax0.axvline(SAMPLING*fLate, ymax = 1., ymin = 0.5, color = f'{colors[ich]}', linestyle = '-')

					ax0.plot(x, BaselineRestoration(wfm = wfm-blines.get("v"), tau = 2900, width = SAMPLING)+blines.get('v')\
						, f'--{colors[ich]}', alpha = 0.8, label = f"rc ch {ich}")
					ax0.axvline(SAMPLING*rcLate, ymax = 0.5, ymin = 0, color = f'{colors[ich]}', linestyle = "--", alpha = 1, label = "RC fixed")
					ax0.axvline(SAMPLING*rcEarly, ymax = 0.5, ymin = 0, color = f'{colors[ich]}',linestyle = "--", alpha = 1)
					
					ax0.axvline(SAMPLING*rcSEarly, ymax = .5, ymin = 0, color = f'{colors[ich]}', linestyle = "dotted", alpha = 1)
					ax0.axvline(SAMPLING*(rcSLate), ymax = 0.5, ymin = 0, color = f'{colors[ich]}', linestyle = "dotted", alpha = 1, label = "RC flat")

		#### done channel loop
		
		#### Event level Processing ####
		tt_evTime[0] = evtime

		#### Fill Tree ####
		tt.Fill()

		##### plotting event level #####
		if args.live == True:
			print("event time %i" % evtime)
			print("blineV", tt_blineV)
			print("blineE", tt_blineE)
			print("fixed charge", tt_qFixed)
			print("flat charge", tt_qFlat)
			print("tt_qRCFixed", tt_qRCFixed)
			print("tt_qRCFlat", tt_qRCFlat)
			print("tt_min", tt_min)
			print("tt_xMin", tt_xMin)

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


	# hFile = gROOT.FindObject(outName)
	# if hFile:
	# 	hFile.Close()

	# with TFile(outName, "RECREATE") as hFile:  
	# 	# h = ROOT.TH1F()  
	# 	# hFile.cd()
	# 	tt.Write()
	# 	print("writing ROOT file to: ", outName, "\n")
	
	hFile = TFile(outName, 'RECREATE')
	#### Closing and Writing #####
	hFile.cd()
	tt.Write()
	hFile.Close()

	print(f"\nDone analysis on {Nevents} Events!\n from run:", fname, "\n")

	print("Wrote output to\n", outName, "\n")

	print(f"% Empty = {100*Nmissing/Nevents/nchans:0.3f}", end = "\n")
	print(f"% Saturated = {100*NSaturated/Nevents/nchans:0.3f}", end = "\n")
	print(f"% Early = {100*Nearly/Nevents/nchans:0.3f}", end = "\n")





