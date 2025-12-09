import os
import sys

CF4HOME = os.environ.get("CF4HOME")
# print("CF4HOME:", CF4HOME)
COMP2PATH = os.path.join(CF4HOME, "comp2")

sys.path.append(os.path.join(CF4HOME, "comp2"))
sys.path.append(os.path.join(CF4HOME, "comp2/pyScripts"))

print("CF4HOME:", CF4HOME)
print("COMP2PATH:", COMP2PATH)