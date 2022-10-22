#!/usr/local/bin/python3

"""
	use this for viewing waveforms. Minimal processing done here. good 'skeleton' for viewing data
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import argparse

import sys
sys.path.append('..')
from SBCcode.DataHandling import ReadBinary as rb

from ozTools import *

print("Austin's makeshift waveform viewer")

# index of recorded waveforms in the event
iSipm = 0 # sipm data is first
iDis = 1 # discriminator data second

if __name__ == "__main__":
	pars = argparse.ArgumentParser()
	pars.add_argument("--first", type = int, required = False, default = 0, help = "what event to start in")
	pars.add_argument("--g", type = int, required = False, default = 0, help = "what pulse group to start on")
	pars.add_argument("--N", type = int, required = False, default = 10, help = "N events to plot together")
	pars.add_argument("--run", type = str, required = False, default = "20210726_0")
	pars.add_argument("--live", type = int, required = False, default = False, help = "do viewing with user input")

	args=pars.parse_args()

	run = args.run # LED data in CF4 gas?
	live = bool(args.live)

	# files and paths
	raw_directory = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/cf4_2022/rawData"

	# setup plotting
	fig, axes = plt.subplots(1, 1, figsize = (8, 6))
	ax0 = axes

	if live == True:
		plt.ion()
		plt.show()

	maxeV = 1 # how many event groups to look through

	# don't loop through pulse groups (badly named event groupings, aka 'events')
	## load the data
	path = os.path.join(raw_directory, run+".bin")
	data = rb.ReadBlock(path)
	
	for key, value in data.items():
		print(key)
		
	# exit()
	# pulse information
	traces = data.get('sipm_traces')
	print(traces.shape)
	print(traces[0][0])
	
	xvals = np.arange(0, len(traces[0][0]) * 16, 16)
	ax0.plot(xvals, traces[0][0])
	plt.show()
	input()
	exit()
	# vScale = data.get('v_scale')
	# offset = data.get('v_offset')

	# loop through 1000 pulses (or less if not many events)
	for ind in range(args.first, args.first + min(len(traces), args.N)):
		
		if ind%1000 == 0: # printing helper
			print("event", ind)

		# recorded sipm waveform
		s1 = traces[ind, iSipm][:-50] # shape
		s1v = s1 * vScale[ind, iSipm] +  offset[ind, iSipm] # shape x scale + offset
		
		d1 = traces[ind, iDis] # discriminator

		# # trigger from discriminator
		# trigger = traces['traces'][ind,1] * traces['v_scale'][ind,1] + traces['v_offset'][ind,1]

		gap = 1 # pointless to multiply by dt, which is just 1ns
		# gap = data.get('dt')[ind, iSipm] * 10**9
		time = np.array(range(len(s1))) * gap #+ data.get("t0")[ind, iSipm]

		ax0.plot(time, s1v, alpha = 0.4)
		# ax0.plot(time, d1)
		# plt.plot(sipm)
		# plt.plot(trigger[:-100]*0.1)

		######## do processing ########
		minTime = np.argmin(s1v)
		blines = getBaselineAvg(data = s1v, start = 0, end = 200)

		if minTime < 300: # pulse is too early
			blines = getBaselineAvg(data = s1v, start = minTime+400, end = minTime+900)
			if live == True:
				print("using end as baseline")

		qRange, rEarly, rLate = minFixRangeIntegrate(data = s1v, bline = blines.get('v'), early = 250, late = 450)
		qSmart1, fEarly, fLate = toFlatRangeIntegrate(data  = s1v, bline = blines.get('v'), blineRMS = blines.get('e'), extra = 0)
			 

		if live == True:
			ax0.clear() # somehow reset but keep last N pulses
			ax0.plot(time, s1v, alpha = 0.4)
			#hold for user input
			print("\n event #", ind)
			print("minTime = ", minTime)
			print("blines = ", blines)
			print("qRange = ", qRange)
			print("qSmart = ", qSmart1)
			plt.pause(0.1)
			plt.show()
			input()
		

	ax0.minorticks_on()
	ax0.grid(True)
	# ax0.set_ylim(-200, 200)

	plt.show()

	# exit()





