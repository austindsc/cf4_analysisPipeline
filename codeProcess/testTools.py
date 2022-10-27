import numpy as np
import matplotlib.pyplot as plt

from ozRCcorrection import *
from ozTools import *

vec = np.arange(0, 10, 1)
vec2 = vec
# a = np.arange(start = 10, end = 100, step = 5)

# print(vec)
# #howAboutThis(vec)
# doesThisCopy(vec)
# print(vec)
# print(vec2)

# print(type(vec))
# print(type(vec2))

t0 = 50.
Vfinal = 5.
Vstart = -1.
trise = 5.
dt = 0.1
tau = 5.

slope = (Vfinal - Vstart)/trise
b = Vstart -slope * (t0 - trise)
line = lambda x : slope * x + b

time = np.arange(0, 1.2*t0, dt)
clean = np.repeat(Vstart, len(time))

for i, t in enumerate(time):

  if t >= t0-trise:
    clean[i] = float(line(t))
  if t >= t0 or line(t) > Vfinal:
    clean[i] = float(Vfinal)
  
time = np.concatenate((time, time+time[-1]))
clean = np.concatenate((clean, np.flip(clean)))

### Done building clean traces ###

rc = np.copy(clean)
for i, vp in enumerate(clean[0:-1]):
  v_i1 = clean[i+1]
  rc[i+1] = rc[i] - (rc[i] - v_i1)*dt/tau


backRC = np.copy(rc)
for i, vp in enumerate(rc[0:-1]):
  
  backRC[i+1] = (rc[i+1] - rc[i]) * tau/dt + rc[i]

###### 
plt.plot(time, clean, '.', label = 'fake data')

plt.plot(time, rc, label = "rc")

# plt.plot(time, 0.5 + backRC, '-', label = "back converted")
# plt.plot(time, 0.5 + undoRC(rc, width = dt, tau = tau), label = "undone")

davidRC = BaselineRestoration(wfm = rc - Vstart, tau = tau, width = dt)
plt.plot(time, davidRC, label = "david's on RC")

plt.legend()

plt.draw()
plt.pause(0.1)
input()
plt.close()
exit()