"""
	A variety of processing tools used to analyze the waveforms. No fitting done here
"""

import numpy as np
import os
import matplotlib.pyplot as plt

from ROOT import ROOT
from ROOT import TH1F

# return a dictionairy with run (date_number) and nEvents
# add on as i require more information
def ozReadFile(path: str, col: int = 0):
	data = np.loadtxt(path, delimiter = '\t', dtype = str, usecols=(col), unpack = True)

	return data

def getBaselineHist(data, vScale):

	high = max(data)
	low = min(data)

	hbline = TH1F( 'hbline', 'baseline histogram', int((high-low)/vScale), low, high)
	for i in data:
		hbline.Fill(i)

	#baseline = hbline.GetBinLowEdge(hbline.GetMaximumBin()) not this one
	rms = hbline.GetRMS()
	mBin = hbline.GetMaximumBin()
	# get average of three bins around peak
	avg = (hbline.GetBinLowEdge(mBin - 1)*hbline.GetBinContent(mBin-1) + hbline.GetBinLowEdge(mBin)*hbline.GetBinContent(mBin) + hbline.GetBinLowEdge(mBin + 1)*hbline.GetBinContent(mBin+1))/(hbline.GetBinContent(mBin-1) + hbline.GetBinContent(mBin) + hbline.GetBinContent(mBin+1))

	return (avg, rms)
	# return avg

def getBaselineAvg(data, start, end):
	return {'v': np.mean(data[start:end]), 'e': np.std(data[start:end])}

# should also include dt in here but whatever
def simpleIntegrate(data, bline):
	#currently setup for negative polarity
	tot = 0
	for i in data:
		height = bline - i
		tot = height + tot
	return tot

def minFixRangeIntegrate(data, bline, early, late):

	iMin = np.argmin(data)
	if iMin - early < 0:
		early = iMin
	if iMin + late > len(data):
		late = len(data) - 1 -iMin

	#currently setup for negative polarity

	tot = 0
	for i in data[iMin - early:iMin + late]:
		height = bline - i
		tot = height + tot
	return tot, iMin-early, iMin + late

# finds the number of array movements until average of +2*nAvg samples is equal to baseline RMS/4 (safety factor)
def findFlat(data, bline, blineRMS, iMin: int, forward: bool = True): # true means go forward
	nAvg = 5

	iFlat = iMin
	avg = np.average(data[iFlat:iFlat+2*nAvg])

	a = 1
	if forward == False:
		a = -1

	while abs(avg - bline) > 0.25*blineRMS and (iFlat > nAvg and iFlat < len(data) - nAvg):
		iFlat = iFlat + a
		avg = np.average(data[iFlat:iFlat+2*nAvg])

	return abs(iFlat-iMin)


# modify this so that it integrates until equal to baseline
# and don't go extra distance on tail end (avoid overshoot on large pulses)
def toFlatRangeIntegrate(data, bline, blineRMS, extra):

	iMin = np.argmin(data)
	forwardGap = findFlat(data = data, bline = bline, blineRMS = blineRMS, iMin = iMin, forward = True)
	reverseGap = findFlat(data = data, bline = bline, blineRMS = blineRMS, iMin = iMin, forward = False)

	early = iMin - reverseGap - extra
	if early < 0:
		early = 0

	late = int(iMin + forwardGap + extra)
	if late > len(data):
		late = len(data) - 1
	# print("flat = ", iMin - reverseGap)
	# print("flat = ", iMin + forwardGap)
	tot = 0
	for i in data[early : late]:
		height = bline - i
		tot = height + tot
	return tot, early, late

def  minAverage(data, bline, spread): # average the data points around the minimum value
	iMin = np.argmin(data)

	return -(np.mean(data[iMin - spread:iMin + spread]) - bline)


