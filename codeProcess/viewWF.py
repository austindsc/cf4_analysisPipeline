#!/usr/local/bin/python3

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

print("Austin's new waveform viewer")

defRun = 202205271113

if __name__ == "__main__":
	pars = argparse.ArgumentParser()
	pars.add_argument("--first", type = int, required = False, default = 0, help = "what event to start in")
	pars.add_argument("--N", type = int, required = False, default = 10, help = "N events to plot together")
	pars.add_argument("--run", type = str, required = False, default = defRun, help = "desired run number")
	pars.add_argument("--live", type = bool, required = False, default = False, help = "do viewing with user input")
	pars.add_argument("--caen", type = bool, required = False, default = True, help = "is this new data (on caen)")

	## for old data type (non CAEN)
	pars.add_argument("--g", type = int, required = False, default = 0, help = "what pulse group to start on. not for CAEN")
	# index of recorded waveforms in the event
	iSipm = 0 # sipm data is first
	iDis = 1 # discriminator data second

	args=pars.parse_args()

	# files and paths
	raw_directory = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/rawData/caen"
	fname = str(args.run)
	
	if not args.caen:
		raw_directory = "old path"
		fname = "other default file and path"

	# setup plotting
	fig, axes = plt.subplots(1, 1, figsize = (8, 6))
	ax0 = axes

	if args.live == True:
		plt.ion()
		plt.show()

	maxeV = 1 # how many event groups to look through

	path = os.path.join(raw_directory, str(args.run)+".bin")
	if not args.caen:
		path = os.path.join(raw_directory, run, str(args.g), "PMTtraces.bin")
	
	## LOAD THE DATA	
	data = rb.ReadBlock(path)
	for key, value in data.items():
		print(key)
	# print(data.get("sample_rate"))

	# OLD data pulse information
	# traces = data.get('traces')
	# vScale = data.get('v_scale')
	# offset = data.get('v_offset')
	
	traces = data.get("sipm_traces") # [event][channel][waveform value]
	times = data.get("time_stamp")
	print(traces.shape)
	print(traces[0][0])

	# loop through 1000 pulses (or less if not many events)
	for ind in range(args.first, args.first + min(len(traces), args.N)):
		
		x0 = np.arange(0, len(traces[ind][0]) * 16, 16)
		x1 = np.arange(0, len(traces[ind][1]) * 16, 16)
		x2 = np.arange(0, len(traces[ind][2]) * 16, 16)
		
		tr0 = traces[ind][0]
		tr1 = traces[ind][1]
		tr2 = traces[ind][2]

		if ind%100 == 0: # printing helper
			print("event", ind)

		# recorded sipm waveform OLD
		# s1 = traces[ind, iSipm][:-50] # shape
		# s1v = s1 * vScale[ind, iSipm] +  offset[ind, iSipm] # shape x scale + offset
		# 
		# d1 = traces[ind, iDis] # discriminator

		# # trigger from discriminator
		# trigger = traces['traces'][ind,1] * traces['v_scale'][ind,1] + traces['v_offset'][ind,1]

		# gap = 1 # pointless to multiply by dt, which is just 1ns
		# gap = data.get('dt')[ind, iSipm] * 10**9
		# time = np.array(range(len(s1))) * gap #+ data.get("t0")[ind, iSipm]
		# OLD
		# ax0.plot(time, d1)
		# plt.plot(sipm)
		# plt.plot(trigger[:-100]*0.1)

		if not args.live:
			ax0.plot(x0, tr0, alpha = 0.4, label = "ch 0")
			ax0.plot(x1, tr1, alpha = 0.4, label = "ch 1")
			ax0.plot(x2, tr2, alpha = 0.4, label = "ch 2")

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

		if args.live == True:
			ax0.clear() # somehow reset but keep last N pulses
			ax0.plot(x0, tr0, alpha = 0.4, label = "ch 0")
			ax0.plot(x1, tr1, alpha = 0.4, label = "ch 1")
			ax0.plot(x2, tr2, alpha = 0.4, label = "ch 2")
			#hold for user input
			# print("\n event #", ind)
			# print("minTime = ", minTime)
			# print("blines = ", blines)
			# print("qRange = ", qRange)
			# print("qSmart = ", qSmart1)
			plt.pause(0.1)
			ax0.minorticks_on()
			ax0.grid(True)
			ax0.legend()
			plt.show()
			input("\n press a key to see next wfm")
		
	if not args.live:
		ax0.minorticks_on()
		ax0.grid(True)
		ax0.legend()
		# ax0.set_ylim(-200, 200)

		plt.show()





