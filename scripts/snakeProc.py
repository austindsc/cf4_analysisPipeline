import subprocess
import sys
import os
import numpy as np

from multiprocessing import Pool

# NRUN = 1000
# defRunList = "/procData/runLists/analysis1/bkgLists/date_0411-48V.txt"
# outDir = "/analysis1_bkg/"

# defRunList = "/procData/runLists/analysis1/bkgLists/date_0414-manyV.txt"
defRunList = "/procData/runLists/analysis2/alphaRuns.txt"
# outDir = "/ana1Co_test/"
outDir = "/analysis2/"

def run_process(command):
	subprocess.run(command, shell = True)
	return

def getMonth(runID):
	return runID[4:4+2]

if __name__ == '__main__':
	pp = str(os.environ["CF4HOME"])

	cA = f"python3 {pp}/codeProcess/newProcVars.py"

	runList = np.loadtxt(pp + defRunList, dtype=str, ndmin = 1)

	###### First check all files exist ##
	allPaths = ["/".join([pp, "rawData", "2022", getMonth(r), r+".bin"]) for r in runList]
	missingPath = []
	for path in allPaths:
		if not os.path.exists(path):
			# print(f"file {path} doesn't exist")
			missingPath.append(path)
	
	if len(missingPath) > 0:
		print("missing paths or files:")
		[print(p) for p in missingPath]
		sys.exit(1)
	else:
		print("all files exist in", defRunList, "\n\n")

	# exit(0)
	processes = []
	iProcs = 4
	# exit(0) # actually run!
	for irun, run in enumerate(runList, start = 1):
		inDir = "/".join(["2022", getMonth(run)])
		
		command = str(cA + " --inDir " + inDir + " --run " + str(run) + " --outDir " + outDir)
		processes.append(command)
		print(command)
		print(irun)

	# running the processes
	pool = Pool(processes = iProcs)
	pool.map(run_process, processes)