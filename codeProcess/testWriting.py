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
"""
print("Austin's new processor")

SATURATED = 6.
SAMPLING = 16.
BLINEEND = 500.
NULLVAL = 0.0
SKIPENDTIME = 500.
defRun = 202205271113
defDir = "2022/02"
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
	savePath = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/"+args.outDir
	if not os.path.exists(savePath):
		print("savePath doesn't exist")
		print(savePath)
		exit()
	fname = str(args.run)
	
	if not args.caen:
		raw_directory = "old path"
		fname = "other default file and path"
	"""
if __name__ == "__main__":
	outName = str("here.root")

	# with TFile("file1.root", "recreate") as outfile:
	hFile = gROOT.FindObject(outName)
	if hFile:
		hFile.Close()

	# with TFile("file1.root", "recreate") as outfile:
	hFile = TFile("file1.root", "RECREATE")
	# 	print("writing ROOT file to: ", outName, "\n")
	
	### defining TTree and all variables ###
	# tt = TTree("tt", "proc static variables, run:")
	# nchans = 3
	h = ROOT.TH1F()
	h.Write()

	# tt_chInfo = np.zeros((nchans), dtype = "int")
	# tt.Branch("chInfo", tt_chInfo, "chInfo[" + str(nchans) +"]/I")
	# pulseDict = {"early":-1, "empty": 0, "normal": 1, "saturated":2	}

	# tt_evTime = np.zeros((1), dtype = "float32")
	# tt.Branch("evTime", tt_evTime, "evTime/F")

	# for i in range(10):
	# 	tt_evTime[0] = i
	# 	tt.Fill()

	# #### Closing and Writing #####
	hFile.cd()
	# tt.Write()
	hFile.Close()






