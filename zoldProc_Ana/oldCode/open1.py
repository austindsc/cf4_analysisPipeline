import numpy as np
import os
import matplotlib.pyplot as plt

import sys
# sys.path.append('..')
from SBCcode.DataHandling import ReadBinary as rb

# 1. load in all runs (date and number)
# 2. get N events? or just try until read fails
# 3. search 

print("hello Austin")

# index of recorded waveforms in the event
iSipm = 0 # sipm data is first
iDis = 1 # discriminator data second

# plt.figure()
fig, axes = plt.subplots(1, 1, figsize = (10, 8))
ax0 = axes
# ax1 = axes[1]

# files and paths
raw_directory = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/atHome/SBC-20-data/"

# normal run (some physics conditions)
run = "20201123_4" # 51V, no saturation, -600mV range
#run = "20201123_0" # 55V, saturated at -200mV

maxeV = 2
maxeV = 0

# loop through pulse groups (badly named events)
for i_eV in range(maxeV + 1):

	## load the data
	path = os.path.join(raw_directory, run, str(i_eV), "PMTtraces.bin")
	data = rb.ReadBlock(path)

	# pulse information
	traces = data.get('traces')
	vScale = data.get('v_scale')
	offset = data.get('v_offset')

	# loop through 1000 pulses
	for ind in range(min(len(traces), 1000)):
		
		if ind%1000 == 0: # printing helper
			print(ind)

		# recorded sipm waveform
		s1 = traces[ind, iSipm] # shape

		s1v = s1 * vScale[ind, iSipm] +  offset[ind, iSipm] # shape x scale + offset

		d1 = traces[ind, iDis] # discriminator

		# # trigger from discriminator
		# trigger = traces['traces'][ind,1] * traces['v_scale'][ind,1] + traces['v_offset'][ind,1]

		gap = 1
		gap = data.get('dt')[ind, iSipm]
		time = np.array(range(len(s1))) * gap + data.get("t0")[ind, iSipm]

		ax0.plot(time, s1v)
		# ax0.plot(time, d1)
		# plt.plot(sipm[:-100])
		# plt.plot(trigger[:-100]*0.1)

ax0.minorticks_on()
ax0.grid(True)
# ax0.set_ylim(-200, 200)

plt.show()

# exit()





