import subprocess
import sys

from multiprocessing import Pool

def run_process(command):
	subprocess.run(command, shell = True)


def getFileList(path):
	

cA = "p311 ../codeProcess/newProcVars.py --N 50000 --outName  "

runList = ["20210723_5", "20210812_2"]

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