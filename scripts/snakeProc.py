import subprocess
import sys

from multiprocessing import Pool

def run_process(command):
	subprocess.run(command, shell = True)

cA = "python3 $CF4HOME/codeProcess/procStaticVars1.py --N 50000 --outName liqMedoutProcStatic2_ "

runList = ["20210726_7", "20210804_4", "20210813_4"]

processes = []
iProcs = 3

for irun, run in enumerate(runList, start = 1):
	command = str(cA + "--run " + str(run))
	processes.append(command)
	print(command)
	print(irun)
	if irun % iProcs == 0:
		pool = Pool(processes = iProcs)
		pool.map(run_process, processes)
		processes = []