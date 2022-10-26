import numpy as np

from ozRCcorrection import *
from ozTools import *

vec = np.arange(0, 10, 1)
vec2 = vec
# a = np.arange(start = 10, end = 100, step = 5)

print(vec)
#howAboutThis(vec)
doesThisCopy(vec)
print(vec)
print(vec2)

print(type(vec))
print(type(vec2))